# Generates a scene containing an elastic cube within a larger cubic enclosure
# Given its elasticity, geometry, and mass, etc. how will the cube bounce?
# Implemented as a deterministic simulation over uncertain properties
# Two ways to provide (fix) observed values - parameters and constraints

# GOAL: Perform sequential inference of elasticity using 3 seconds of observations, followed by 3 seconds of prediction
# GOAL: Test what simulation settings (i.e. resource-rational approximations) best match human performance:
#        * using a rigid body
#        * using a sphere
#        * adding noise to orientation
# DONE: Print / write the traces to a .csv file, check whether it recovers elasticity
# TODO: Set initial position using RealFlow data
# TODO: Visualize by plotting bounce locations in 3D, with walls 'sketched' in
# TODO: Add ceiling and fourth wall; consider different restitution values
# TODO: Multiple forward passes per particle, with some noise added over velocity

using Accessors
using Distributions
using Gen
using PyCall
using PhySMC
using PhyBullet
using DataFrames
using CSV

bullet = pyimport("pybullet")
pybullet_data = pyimport("pybullet_data")
include("truncatednorm.jl")
include("plots.jl")


# Generative Model

# Samples initial latent estimates
@gen function sample_latents(latents::RigidBodyLatents)
    mass = {:mass} ~ gamma(1.2, 10.)
    res = {:restitution} ~ uniform(0, 1)
    new_latents = RigidBodyLatents(setproperties(latents.data, mass=mass, restitution=res))
    return new_latents
end

# Sets initial scene configuration
# In the future this will be an inverse graphics module
@gen function generate_scene(sim::PhySim, init_position::Vector{Float64}, init_velocity::Vector{Float64}, mass::Float64=1.0, restitution::Float64=0.9)

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
    # Sample initial position
    init_position_x = {:init_state => :x0} ~ uniform(-1, 1)
    init_position_y = {:init_state => :y0} ~ uniform(-1, 1)
    init_position_z = {:init_state => :z0} ~ uniform(1, 3)
    start_position = [init_position_x, init_position_y, init_position_z]
    =#

    startPosition = init_position
    startOrientation = bullet.getQuaternionFromEuler([0, 0, 1])

    # Create and position cube
    cubeBody = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.2, .2, .2])
    cube = bullet.createMultiBody(baseCollisionShapeIndex=cubeBody, basePosition=startPosition, baseOrientation=startOrientation)

    #=
    # Create and position sphere
    sphereBody = bullet.createCollisionShape(bullet.GEOM_SPHERE, radius=.2)
    sphere = bullet.createMultiBody(baseCollisionShapeIndex=sphereBody, basePosition=startPositionCube, baseOrientation=startOrientationCube)
    =#

    # Set kinematics and dynamics dynamics
    bullet.changeDynamics(cube, -1, mass=mass, restitution=restitution)
    bullet.resetBaseVelocity(cube, linearVelocity=init_velocity)

    # Store representation of cube's initial state
    init_state = BulletState(sim, [RigidBody(cube)])

    return init_state
end

# Adds measurement noise to ground truth position
@gen function generate_observation(k::RigidBodyState)
    obs = {:position} ~ broadcasted_normal(k.position, 0.1)
    return obs
end

# Syncs Bullet to the given state, generates the next state, then samples an observation
@gen function kernel(t::Int, current_state::BulletState, sim::BulletSim)

    # Synchronizes state, then use Bullet to generate next state
    next_state::BulletState = PhySMC.step(sim, current_state)

    # Applies noise to x, y, and z positions
    {:observation} ~ Gen.Map(generate_observation)(next_state.kinematics)

    return next_state
end

# Samples latents from their priors, then runs complete forward simulation
@gen function generate_trajectory(T::Int, init_state::BulletState, sim::BulletSim)

    # Resample the target object's latents - mass and restitution
    latents = {:latents} ~ Gen.Map(sample_latents)(init_state.latents)
    
    # init_state = generate_scene(sim, latents.mass, latents.restitution)
    init_state = Accessors.setproperties(init_state; latents=latents)

    # Simulate T time steps
    states = {:trajectory} ~ Gen.Unfold(kernel)(T, init_state, sim)

    return states
end


# Inference

# Parses the full sequence of observations from choice map
function get_observations(choices::Gen.ChoiceMap, T::Int)
    observations = Vector{Gen.ChoiceMap}(undef, T)
    for i=1:T
        cm = choicemap()
        addr = :trajectory => i => :observation => :position
        set_submap!(cm, addr, get_submap(choices, addr))
        observations[i] = cm
    end
    return observations
end

# Propose new latents - for MCMC
@gen function proposal(trace::Gen.Trace)

    # Read current mass and restitution estimates from trace
    choices = get_choices(trace)
    prev_mass = choices[:latents => 1 => :mass]
    prev_res  = choices[:latents => 1 => :restitution]

    # Resample mass and restitution from Gaussians centered at their previous values
    mass = {:latents => 1 => :mass} ~ trunc_norm(prev_mass, 0.1, 0.0, Inf)
    restitution = {:latents => 1 => :restitution} ~ trunc_norm(prev_res, 0.1, 0.0, 1.0)

    return (mass, restitution)
end

# Particle filter
function infer(gm_args::Tuple, obs::Vector{Gen.ChoiceMap}, num_particles::Int=20)

    # Update first argument (time step) to generative model
    model_args(t) = (t, gm_args[2:3]...)

    # Initiliaze the particle filter with no observations
    state = Gen.initialize_particle_filter(generate_trajectory, model_args(0), EmptyChoiceMap(), num_particles)

    # Iteratively simulate and filter particles
    argdiffs = (UnknownChange(), NoChange(), NoChange())
    for (t, obs) = enumerate(obs)
        @elapsed begin
            
            # Decide whether to cull and resample poorly performing particles
            Gen.maybe_resample!(state, ess_threshold=num_particles/2)

            # Simulate the next time step and score against new observation
            Gen.particle_filter_step!(state, model_args(t), argdiffs, obs)

            # Rejuvenation
            for i=1:num_particles
                state.traces[i], _ = Gen.mh(state.traces[i], proposal, ())
            end

        end
    end

    return Gen.sample_unweighted_traces(state, 5)
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

function write_to_csv(particles, fname=joinpath(pwd(), "test.csv"))

    println("Writing simulation data to " * fname)
    particle_data = DataFrame(particle=Int[], elasticity=[], frame=Int[], x=[], y=[], z=[], ox=[], oy=[], oz=[], ow=[])

    for (p, particle) in enumerate(particles)
        ela = particle[:latents => 1 => :restitution]
        for (f, frame) in enumerate(particle[:trajectory])
            pos = convert(Vector, frame.kinematics[1].position)
            ori = convert(Vector, frame.kinematics[1].orientation)
            data = [p; ela; f; pos; ori]
            push!(particle_data, data)
        end
    end

    CSV.write(fname, particle_data)
end


# Tests

# Tests whether a range of elasticity settings produce plausible trajectories
function test_elasticity()
    args = (60, sim, init_state)
    gt_constraints = choicemap((:latents => 1 => :restitution, 0.8), (:latents => 1 => :mass, 1.0))
    trace = first(generate(simulation, args, gt_constraints))

    traces = [trace]
    for i=2:5
        res = i*0.2
        gt_constraints = choicemap((:latents => 1 => :restitution, res), (:latents => 1 => :mass, 1.0))
        trace = first(generate(simulation, args, gt_constraints))
        push!(traces, trace)
    end

    for trace in traces
        choices = get_choices(trace)
        display(choices)
        positions = [choices[:trajectory => i => :observation => 1 => :position] for i=35:60]
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

    for i=1:5
        # Read ground truth initial velocity
        fname = "RealFlowData/Cube_Ela9_Var118.csv"
        data = CSV.read(fname, DataFrame)
        datum = values(data[1, 11:13])
        initial_velocity = [datum[1], datum[3], datum[2]]
        println(initial_velocity)

        # Read ground truth initial position
        fname = "RealFlowData/Cube_Ela9_Var118_observed.csv"
        data = CSV.read(fname, DataFrame)
        datum = values(data[1, :])
        initial_position = [datum[1], datum[3], datum[2]]
        println(initial_position)

        init_state = generate_scene(sim, initial_position, initial_velocity)

        # Read ground truth trajectory
        observations = Vector{Gen.ChoiceMap}(undef, size(data)[1]-1)
        zs = Vector{Float64}(undef, size(data)[1]-1)
        for i=1:size(data)[1]-1
            addr = :trajectory => i => :observation => :position
            datum = values(data[i, :])
            new_datum = [datum[1], datum[3], datum[2]]
            cm = Gen.choicemap((addr, new_datum))
            observations[i] = cm
            zs[i] = datum[2]
        end
        gif(animate_observations(zs))

        #=
        # Generate ground truth trajectory
        gt_constraints = choicemap((:latents => 1 => :restitution, 0.7), (:latents => 1 => :mass, 1.0))
        ground_truth = first(generate(generate_trajectory, args, gt_constraints))

        # Display ground truth trajectory
        gif(animate_trace(ground_truth), fps=24)
        gt_choices = get_choices(ground_truth)
        display(gt_choices)
        observations = get_observations(gt_choices, args[1])
        =#

        # Infer elasticity from observed trajectory 
        args = (29, init_state, sim)
        result = infer(args, observations)
        display(get_choices(result[1]))
        gif(animate_traces(result), fps=24)
        write_to_csv(result, "BulletData/observations.csv")
        
        # For each particle, predict the next 60 time steps
        ppd = predict(result, 90)    
        gif(animate_traces(ppd), fps=24)  
        write_to_csv(ppd, "BulletData/predictions.csv")   
    end

    bullet.disconnect()
end


main()
