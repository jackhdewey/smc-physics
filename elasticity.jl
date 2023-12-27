# Generates a scene containing an elastic cube within a larger cubic enclosure
# Given its elasticity, geometry, and mass, etc. how will the cube bounce?
# Implemented as a deterministic simulation over uncertain properties
# Two ways to provide (fix) observed values - parameters and constraints

# GOAL: Perform sequential inference of elasticity using 3 seconds of observations, followed by 3 seconds of prediction
# GOAL: Test what simulation settings (i.e. resource-rational approximations) best match human performance:
#        * using a rigid body
#        * using a sphere
#        * adding noise to orientation
# TODO: Visualize by plotting bounce locations in 3D, with walls 'sketched' in
# TODO: Add ceiling and fourth wall; consider different restitution values
# TODO: Multiple forward passes per particle, with some noise added over velocity
# DONE: Print / write the traces to a .csv file, check whether it recovers elasticity
# DONE: Add a prediction phase that simulates forward some output particles
# DONE: Try out elasticity values 0.2 - 1.0
# DONE: Try a sphere

using Accessors
using Distributions
using Gen
using PyCall
using PhySMC
using PhyBullet
using DataFrames
using CSV

using DataFrames
using CSV

bullet = pyimport("pybullet")
pybullet_data = pyimport("pybullet_data")
include("truncatednorm.jl")
include("plots.jl")


# Generative Model

# Updates latent variables
# The latents are the unobserved variables
# They need to be updated in the mental representation
function update_latents(latents::RigidBodyLatents, mass::Float64, res::Float64)
    RigidBodyLatents(setproperties(latents.data, mass=mass, restitution=res))
end

# Samples initial latent estimates
@gen function sample_from_prior(latents::RigidBodyLatents)
    mass = {:mass} ~ gamma(1.2, 10.)
    res = {:restitution} ~ uniform(0, 1)
    new_latents = update_latents(latents, mass, res)
    return new_latents
end

# Adds measurement noise to ground truth position
@gen function generate_observation(k::RigidBodyState)
    obs = {:position} ~ broadcasted_normal(k.position, 0.1)
    return obs
end

# Sets initial scene configuration
# In the future this will be an inverse graphics module
function generate_scene(sim::PhySim, mass::Float64=1.0, restitution::Float64=0.9)

    bullet.setGravity(0, 0, -10)

    # Create and position ground plane
    planeID = bullet.createCollisionShape(bullet.GEOM_PLANE)
    plane = bullet.createMultiBody(baseCollisionShapeIndex=planeID, basePosition=[0, 0, 0])
    bullet.changeDynamics(plane, -1, mass=0.0, restitution=0.5)

    # Create and position walls
    wall1ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.02, 1, 1])
    quaternion = bullet.getQuaternionFromEuler([0, 0, 0])
    plane2 = bullet.createMultiBody(baseCollisionShapeIndex=wall1ID, basePosition=[-1, 0, 1], baseOrientation=quaternion)
    bullet.changeDynamics(plane2, -1, mass=0.0, restitution=0.5)

    wall2ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[1, .02, 1])
    plane3 = bullet.createMultiBody(baseCollisionShapeIndex=wall2ID, basePosition=[0, 1, 1], baseOrientation=quaternion)
    bullet.changeDynamics(plane3, -1, mass=0.0, restitution=0.5)

    wall3ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.02, 1, 1])
    plane4 = bullet.createMultiBody(baseCollisionShapeIndex=wall3ID, basePosition=[1, 0, 1], baseOrientation=quaternion)
    bullet.changeDynamics(plane4, -1, mass=0.0, restitution=0.5)

    #=
    # Create and position cube
    cubeBody = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.2, .2, .2])
    startOrientationCube = bullet.getQuaternionFromEuler([0, 0, 1])
    cube = bullet.createMultiBody(baseCollisionShapeIndex=cubeBody, basePosition=[0., 0., 3.], baseOrientation=startOrientationCube)
    bullet.changeDynamics(cube, -1, mass=mass, restitution=restitution)
    =#

    # Create and position sphere
    sphereBody = bullet.createCollisionShape(bullet.GEOM_SPHERE, radius=.2)
    startOrientationCube = bullet.getQuaternionFromEuler([0, 0, 1])
    sphere = bullet.createMultiBody(baseCollisionShapeIndex=sphereBody, basePosition=[0., 0., 3.], baseOrientation=startOrientationCube)
    bullet.changeDynamics(sphere, -1, mass=mass, restitution=restitution)

    # Generate representation of sphere's initial state
    rb_sphere = RigidBody(sphere)
    init_state = BulletState(sim, [rb_sphere])

    return init_state
end

# Deterministically generates the next state, then samples an observation
@gen function kernel(t::Int, current_state::BulletState, sim::BulletSim)

    # Synchronizes state and asks Bullet to generate next time step
    next_state::BulletState = PhySMC.step(sim, current_state)

    # Applies noise independently to x, y, and z position
    {:observation} ~ Gen.Map(generate_observation)(next_state.kinematics)

    return next_state
end

# Samples latents from priors, then run complete forward simulation
@gen function generate_scene_trajectory(T::Int, sim::BulletSim, init_state::BulletState)

    latents = {:latents} ~ Gen.Map(sample_from_prior)(init_state.latents)
    init_state = Accessors.setproperties(init_state; latents=latents)

    states = {:trajectory} ~ Gen.Unfold(kernel)(T, init_state, sim)

    return states
end


# Inference

function get_observations(trace, T::Int)
    observations = Vector{Gen.ChoiceMap}(undef, T)
    for i=1:T
        obs = choicemap()
        address = :kernel => i => :observation
        set_submap!(obs, address, get_submap(gt_choices, address))
        observations[i] = obs
    end
end

# Propose new latents
@gen function proposal(trace::Gen.Trace)

    # Read current mass and restitution estimates from trace
    choices = get_choices(trace)
    prev_mass = choices[:prior => 1 => :mass]
    prev_res  = choices[:prior => 1 => :restitution]

    # Resample mass and restitution from Gaussians centered at their previous values
    mass = {:prior => 1 => :mass} ~ trunc_norm(prev_mass, 0.1, 0.0, Inf)
    restitution = {:prior => 1 => :restitution} ~ trunc_norm(prev_res, 0.1, 0.0, 1.0)

    return (mass, restitution)
end

# Particle filter
function infer(gm_args::Tuple, obs::Vector{Gen.ChoiceMap}, num_particles::Int=20)

    # Function to update first argument (time step) to generative model
    model_args(t) = (t, gm_args[2:3]...)

    # Initiliaze the particle filter with no observations
    state = Gen.initialize_particle_filter(simulation, model_args(0), EmptyChoiceMap(), num_particles)
    argdiffs = (UnknownChange(), NoChange(), NoChange())

    for (t, obs) = enumerate(obs)
        @elapsed begin

            # Rejuvenation
            for i=1:num_particles
                state.traces[i], _ = Gen.mh(state.traces[i], proposal, ())
            end
            
            # Decide whether to cull poorly performing particles
            Gen.maybe_resample!(state, ess_threshold=num_particles/2)

            # Simulate the next time step and score against observation
            Gen.particle_filter_step!(state, model_args(t), argdiffs, obs)
        end
    end

    return state.traces
end


# Prediction

function predict(result, T::Int) 

    n = length(result)
    ppd = Vector{Gen.Trace}(undef, n)
    for i=1:n
        prev_args = get_args(result[i])
        new_args = (prev_args[1] + T, prev_args[2:end]...)
        arg_diffs = (UnknownChange(), NoChange(), NoChange())
        ppd[i] = first(Gen.update(result[i], new_args, arg_diffs, choicemap()))
    end

    #=
    for trace in result
        args = (60, sim, BulletState(sim, [rb_cube]))
        constraints = Gen.choicemap()
        constraints[:prior => 1 => :mass] = trace[:prior => 1 => :mass]
        constraints[:prior => 1 => :restitution] = trace[:prior => 1 => :restitution]
        Gen.generate(simulation, args, constraints)
    end
    =#

    return ppd
end


# Data

function write_to_csv(coll, fname=joinpath(pwd(), "test.csv"))

    println("Writing simulation data to " * fname)
    particles_data = DataFrame(particle=Int[], frame=Int[], x=[], y=[], z=[], ox=[], oy=[], oz=[], ow=[])
    fname = "test.csv"

    for p_num = 1:length(coll)
        particle = coll[p_num]
        for (f, frame) in enumerate(particle[:kernel])
            pos = convert(Vector, frame.kinematics[1].position)
            ori = convert(Vector, frame.kinematics[1].orientation)
            data = [p_num; f; pos; ori]

            push!(particles_data, data)
        end
    end

    CSV.write(fname, particles_data)
end

function read_from_csv(fname)
    data = CSV.read(fname, DataFrame)
end


# Tests

# Tests whether a range of elasticity settings produce plausible trajectories
function test_elasticity()
    args = (60, sim, init_state)
    gt_constraints = choicemap((9:prior => 1 => :restitution, 0.8), (:prior => 1 => :mass, 1.0))
    trace = first(generate(simulation, args, gt_constraints))

    traces = [trace]
    for i=2:5
        res = i*0.2
        gt_constraints = choicemap((:prior => 1 => :restitution, res), (:prior => 1 => :mass, 1.0))
        trace = first(generate(simulation, args, gt_constraints))
        push!(traces, trace)
    end

    for trace in traces
        choices = get_choices(trace)
        display(choices)
        positions = [choices[:kernel => i => :observation => 1 => :position] for i=35:60]
        zs = map(x -> x[3], positions)
        z_max = maximum(zs)
        coefficient_of_restitution = z_max / 3.0
        println(coefficient_of_restitution)
    end
end


# Main

function main()

    # Initialize simulation context 
    client = bullet.connect(bullet.DIRECT)::Int64
    bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
    #bullet.resetSimulation(bullet.RESET_USE_DEFORMABLE_WORLD)
    bullet.resetDebugVisualizerCamera(5, 0, -5, [0, 0, 2])
    sim = BulletSim(step_dur=1/30; client=client)

    # Generate ground truth trajectory
    init_state = generate_scene(sim)
    args = (60, sim, init_state)
    gt_constraints = choicemap((:latents => 1 => :restitution, 0.7), (:latents => 1 => :mass, 1.0))
    ground_truth = first(generate(generate_scene_trajectory, args, gt_constraints))

    gif(animate_trace(ground_truth), fps=24)
    display(get_choices(ground_truth))

    # TODO: Add option to read ground truth from a .csv 

    # Infer elasticity from position observations
    get_observations(ground_truth, args[1])
    result = infer(args, observations)
    write_to_csv(result)
    gif(animate_traces(result), fps=24)
    
    # For each particle, predict the next 60 time steps
    ppd = predict(result, 60)    
    gif(animate_traces(ppd), fps=24)        

    bullet.disconnect()
end


main()
