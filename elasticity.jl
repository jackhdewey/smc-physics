# Generates a scene containing an elastic cube within a larger cubic enclosure
# Given its elasticity, geometry, and mass, etc. how will the cube bounce?
# Implemented as a deterministic simulation over uncertain properties
# Two ways to provide (fix) observed values - parameters and constraints

# GOAL: Perform sequential inference of elasticity using 3 seconds of observations, followed by 3 seconds of prediction
# GOAL: Test what simulation settings (i.e. resource-rational approximations) best match human performance:
#        * using a rigid body
#        * using a sphere
#        * adding noise to orientation
# TODO: Add ceiling and fourth wall; consider different restitution values
# TODO: Print / write the traces to a .csv file, check whether it recovers elasticity
# DONE: Try out elasticity values 0.2 - 1.0
# DONE: Try a sphere

using Accessors
using Distributions
using Gen
using PyCall
using PhySMC
using PhyBullet
bullet = pyimport("pybullet")
pybullet_data = pyimport("pybullet_data")
include("truncatednorm.jl")
include("plots.jl")


# Utilities

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

# Rejuvenate latent estimates
# For using during MCMC
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

function write_to_csv(coll, fname=joinpath(pwd(), "test.csv"))
    particles_data = []
    append!(particles_data, ["particle, timestep, x, y, z, qx, qy, qz, qw"])

    for p_num = 1:length(coll)
        particle = coll[p_num]
        for (t, frame) in enumerate(particle[:kernel])

            pos = convert(Vector, frame.kinematics[1].position)
            ori = convert(Vector, frame.kinematics[1].orientation)

            data = Any[p_num; t; pos; ori]

            push!(particles_data, data)
        end
    end

    open(fname, "w") do f
        println(f, particles_data[1])
        for row in particles_data[2:end]
            println(f, join(row, ","))
        end
    end
end

# Generative Model

# Sets initial scene configuration
# In the future this will be an inverse graphics module
function init_scene(mass::Float64=1.0, restitution::Float64=0.9)

    bullet.setGravity(0, 0, -10)

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

    return sphere
end

# Adds measurement noise to ground truth position
@gen function generate_observation(k::RigidBodyState)
    pos = k.position
    obs = {:position} ~ broadcasted_normal(pos, 0.1)
    return obs
end

# Deterministically generates the next state, then samples an observation
@gen function kernel(t::Int, current_state::BulletState, sim::BulletSim)
    next_state::BulletState = PhySMC.step(sim, current_state)
    {:observation} ~ Gen.Map(generate_observation)(next_state.kinematics)
    return next_state
end

# Samples latents from priors, then runs complete forward simulation
@gen function simulation(T::Int, sim::BulletSim, init_state::BulletState)
    latents = {:prior} ~ Gen.Map(sample_from_prior)(init_state.latents)
    init_state = Accessors.setproperties(init_state; latents=latents)
    states = {:kernel} ~ Gen.Unfold(kernel)(T, init_state, sim)
    return states
end


# Inference

# Particle filter
function infer(gm_args::Tuple, obs::Vector{Gen.ChoiceMap}, num_particles::Int=20)

    # Function to update first argument (time step) to generative model
    model_args(t) = (t, gm_args[2:3]...)

    # Initiliaze the particle filter with no observations
    state = Gen.initialize_particle_filter(simulation, model_args(0), EmptyChoiceMap(), num_particles)
    argdiffs = (UnknownChange(), NoChange(), NoChange())

    for (t, obs) = enumerate(obs)
        @elapsed begin
            for i=1:num_particles
                state.traces[i], _ = Gen.mh(state.traces[i], proposal, ())
            end
            
            Gen.maybe_resample!(state, ess_threshold=num_particles/2)

            Gen.particle_filter_step!(state, model_args(t), argdiffs, obs)
        end
    end

    # TODO: Add a prediction phase that simulates forward some output particles
    # TODO: Multiple forward passes per particle, with some noise added over velocity

    return state.traces
end

# Prediction

function predict() 

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
    client = bullet.connect(bullet.DIRECT)::Int64
    bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
    bullet.resetDebugVisualizerCamera(5, 0, -5, [0, 0, 2])
    # bullet.resetSimulation(bullet.RESET_USE_DEFORMABLE_WORLD)

    # Initialize simulation context 
    sim = BulletSim(; client=client)

    # Initialize scene
    rb_cube = RigidBody(init_scene())
    init_state = BulletState(sim, [rb_cube])

    # Generate ground truth trajectory
    args = (60, sim, init_state)
    gt_constraints = choicemap((:prior => 1 => :restitution, 0.7), (:prior => 1 => :mass, 1.0))
    ground_truth = first(generate(simulation, args, gt_constraints))
    gif(animate_trace(ground_truth), fps=24)

    # TODO: Add option to read ground truth from a .csv
   
    # Transfer position observations from trace to a vector
    gt_choices = get_choices(ground_truth)
    display(gt_choices)

    t = args[1]
    observations = Vector{Gen.ChoiceMap}(undef, t)
    for i = 1:t
        obs = choicemap()
        address = :kernel => i => :observation
        set_submap!(obs, address, get_submap(gt_choices, address))
        observations[i] = obs
    end

    # Infer ground truth elasticity
    result = infer(args, observations)
    display(get_choices(result[1]))
    println(length(result))

    #=
    for trace in result
        args = (60, sim, BulletState(sim, [rb_cube]))
        constraints = Gen.choicemap()
        constraints[:prior => 1 => :mass] = trace[:prior => 1 => :mass]
        constraints[:prior => 1 => :restitution] = trace[:prior => 1 => :restitution]
        Gen.generate(simulation, args, constraints)
    end
    =#

    # TODO: Write elasticities and trajectories for each particle to .csv 

    # Visualize particles
    gif(animate_traces(result), fps=24)
    write_to_csv(result)

    bullet.disconnect()
end


main()
