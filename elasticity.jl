# Generates a scene containing an elastic cube or sphere within a larger cubic enclosure
# Given its elasticity and geometry, how will the object bounce (what will be its trajectory)?
# Implemented as a deterministic simulation over uncertain properties

# GOAL: Perform sequential inference of elasticity using 3 seconds of observations, followed by 3 seconds of prediction
# GOAL: Test what simulation settings (i.e. resource-rational approximations) best match human performance:
#           * using a rigid body
#           * using a sphere
#           * adding noise to orientation
#           * etc.
#
# TODO: Fix off-by-one error
#
# DONE: Set initial position and velocity using RealFlow data
# DONE: Print / write the traces to a .csv file, check whether it recovers elasticity
# DONE: Visualize by plotting bounce locations in 3D
# DONE: Save itermediate particle filter states
#
# CONSIDER: Changing prior for elasticity
# CONSIDER: Adding noise over velocity
# CONSIDER: Handling collisions separately

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
include("Utilities/truncatednorm.jl")
include("Utilities/plots.jl")
include("Utilities/data.jl")


# Generative Model

# Sets initial scene configuration, including kinematics and dynamics
# In the future this will be an inverse graphics module
function generate_scene(sim::PhySim, init_position::Vector{Float64}, init_velocity::Vector{Float64}, mass::Float64=1.0, restitution::Float64=0.9)

    # Create and position ground plane
    planeID = bullet.createCollisionShape(bullet.GEOM_PLANE)
    plane = bullet.createMultiBody(baseCollisionShapeIndex=planeID, basePosition=[0, 0, 0])
    bullet.changeDynamics(plane, -1, mass=0.0, restitution=0.5)

    # Create and position walls
    wall1ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.01, 0.5, 0.5])
    quaternion = bullet.getQuaternionFromEuler([0, 0, 0])
    wall1 = bullet.createMultiBody(baseCollisionShapeIndex=wall1ID, basePosition=[-0.5, 0, 0.5], baseOrientation=quaternion)
    bullet.changeDynamics(wall1, -1, mass=0.0, restitution=0.5)

    wall2ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.5, 0.01, 0.5])
    wall2 = bullet.createMultiBody(baseCollisionShapeIndex=wall2ID, basePosition=[0, 0.5, 0.5], baseOrientation=quaternion)
    bullet.changeDynamics(wall2, -1, mass=0.0, restitution=0.5)

    wall3ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.01, 0.5, 0.5])
    wall3 = bullet.createMultiBody(baseCollisionShapeIndex=wall3ID, basePosition=[0.5, 0, 0.5], baseOrientation=quaternion)
    bullet.changeDynamics(wall3, -1, mass=0.0, restitution=0.5)

    wall4ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.5, 0.01, 0.5])
    wall4 = bullet.createMultiBody(baseCollisionShapeIndex=wall4ID, basePosition=[0, -0.5, 0.5], baseOrientation=quaternion)
    bullet.changeDynamics(wall4, -1, mass=0.0, restitution=0.5)
    
    ceilingID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.5, 0.5, 0.01])
    ceiling = bullet.createMultiBody(baseCollisionShapeIndex=ceilingID, basePosition=[0, 0, 1], baseOrientation=quaternion)
    bullet.changeDynamics(ceiling, -1, mass=0.0, restitution=0.5)

    # Set initial position and orientation
    startPosition = init_position
    startOrientation = bullet.getQuaternionFromEuler([0, 0, 1])

    #=
    # Create, position, and orient cube
    cubeBody = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.05, .05, .05])
    cube = bullet.createMultiBody(baseCollisionShapeIndex=cubeBody, basePosition=startPosition, baseOrientation=startOrientation)
    =#

    # Create, position, and orient sphere
    sphereBody = bullet.createCollisionShape(bullet.GEOM_SPHERE, radius=.1)
    sphere = bullet.createMultiBody(baseCollisionShapeIndex=sphereBody, basePosition=startPosition, baseOrientation=startOrientation)

    # Set dynamic properties and initial kinematic state
    bullet.setGravity(0, 0, -10)
    bullet.changeDynamics(sphere, -1, mass=mass, restitution=restitution)
    bullet.resetBaseVelocity(sphere, linearVelocity=init_velocity)

    # Store representation of cube's initial state
    init_state = BulletState(sim, [RigidBody(sphere)])

    return init_state
end

# Adds measurement noise to ground truth position
@gen function generate_observation(k::RigidBodyState)
    obs = {:position} ~ broadcasted_normal(k.position, 0.1)
    return obs
end

# Samples an observation, syncs Bullet to the given state, and generates the next state
@gen function kernel(t::Int, current_state::BulletState, sim::BulletSim)

    # Applies noise to x, y, and z positions
    {:observation} ~ Gen.Map(generate_observation)(current_state.kinematics)

    # Synchronizes state, then use Bullet to generate next state
    next_state::BulletState = PhySMC.step(sim, current_state)

    return next_state
end

# Samples initial latent estimate
@gen function sample_latents(latents::RigidBodyLatents)
    # mass = {:mass} ~ gamma(1.2, 10.)
    res = {:restitution} ~ uniform(0, 1)
    new_latents = RigidBodyLatents(setproperties(latents.data, restitution=res))
    return new_latents
end

# Given an initial state, samples latents from their priors then runs a complete forward simulation
@gen function generate_trajectory(T::Int, init_state::BulletState, sim::BulletSim)

    # Sample values for the target object's latents - i.e. a restitution
    latents = {:latents} ~ Gen.Map(sample_latents)(init_state.latents)
    init_state = Accessors.setproperties(init_state; latents=latents)

    # Simulate T time steps
    states = {:trajectory} ~ Gen.Unfold(kernel)(T, init_state, sim)

    return states
end

# Generates a single trace with the specified elasticity and arguments
function generate_ground_truth(args::Tuple, restitution::Float64)
  
    gt_constraints = choicemap((:latents => 1 => :restitution, restitution))
    ground_truth = first(generate(generate_trajectory, args, gt_constraints))
        
    # Display trajectory
    gt_choices = get_choices(ground_truth)
    display(gt_choices)
    gif(animate_trace(ground_truth), fps=24)

    return ground_truth
end

# Parses a sequence of observations from a choice map
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


# Inference and Prediction

# Resamples latent(s) from a Gaussian proposal distribution
@gen function proposal(trace::Gen.Trace)

    choices = get_choices(trace)

    # Resample restitution from a Gaussian centered at the previous estimate
    prev_res  = choices[:latents => 1 => :restitution]    
    restitution = {:latents => 1 => :restitution} ~ trunc_norm(prev_res, 0.1, 0.0, 1.0)

    # prev_mass = choices[:latents => 1 => :mass]
    # mass = {:latents => 1 => :mass} ~ trunc_norm(prev_mass, 0.1, 0.0, Inf)

    return (restitution)
end

# Particle filter
function infer(fname, gm_args::Tuple, obs::Vector{Gen.ChoiceMap}, num_particles::Int=20)

    # Extract test information    
    tokens = split(fname, "_")

    # Initiliaze the particle filter (with no observations)
    state = Gen.initialize_particle_filter(generate_trajectory, (0, gm_args[2:3]...), EmptyChoiceMap(), num_particles)
    
    # Iteratively simulate and filter particles
    argdiffs = (UnknownChange(), NoChange(), NoChange())
    for (t, obs) = enumerate(obs)
        @elapsed begin
            
            # Decide whether to cull and resample poorly performing particles
            Gen.maybe_resample!(state, ess_threshold=num_particles/2)

            # Simulate the next time step and score against new observation
            Gen.particle_filter_step!(state, (t, gm_args[2:3]...), argdiffs, obs)

            # Resample elasticity, test by simulating to current time step, then accept / reject
            for i=1:num_particles
                state.traces[i], _ = Gen.mh(state.traces[i], proposal, ())
            end
            
            # Dump current particles to a .csv file
            fname = string("BulletData/Intermediate/particles_", t, "_", tokens[2], "_", tokens[3])
            write_to_csv(state.traces, fname)

        end
    end

    # Gen.sample_unweighted_traces(state, 5)
    weights = Gen.get_log_weights(state)

    return state.traces, weights
end

# Given the inferred elasticity, predict the next T observations
function predict(particles, T::Int) 

    # For each particle
    n = length(particles)
    ppd = Vector{Gen.Trace}(undef, n)
    for i=1:n
        # Reset the number of time steps, leave other arguments and existing trace intact
        prev_args = get_args(particles[i])
        new_args = (prev_args[1] + T, prev_args[2:end]...)
        arg_diffs = (UnknownChange(), NoChange(), NoChange())
        
        # Generate the next T time steps
        ppd[i] = first(Gen.update(particles[i], new_args, arg_diffs, choicemap()))
    end

    return ppd
end


# Main

function main()

    # Read ground truth trajectories
    dir = "RealFlowData/"
    fnames = readdir(dir)
    fnames = filter_unwanted_filenames(fnames)
    println(fnames)

    for i in eachindex(fnames)

        # Initialize simulation context 
        client = bullet.connect(bullet.DIRECT)::Int64
        bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
        #bullet.resetSimulation(bullet.RESET_USE_DEFORMABLE_WORLD)
        bullet.resetDebugVisualizerCamera(4, 0, -4, [0, 0, 2])
        sim = BulletSim(step_dur=1/30; client=client)

        # Initialize simulation using observed data
        fname = fnames[i]
        initial_position, initial_velocity, observations = read_observation_file(fname)
        init_state = generate_scene(sim, initial_position, initial_velocity)
        args = (30, init_state, sim)

        # Infer elasticity from observed trajectory 
        results, weights = infer(fname, args, observations, 20)
        gif(animate_traces(results), fps=24)

        # Write inferred elasticities to a .csv file
        tokens = split(fname, "_")
        fname = string("BulletData/Observations/observations_", tokens[2], "_", tokens[3])
        write_to_csv(results, fname)
            
        # For each particle, predict the next 90 time steps
        ppd = predict(results, 90)
        gif(animate_traces(ppd), fps=24)

        # Write predicted trajectory to a .csv file
        fname = string("BulletData/Predictions/predictions_", tokens[2], "_", tokens[3])
        write_to_csv(ppd, fname)

        bullet.disconnect()
    end
end


main()
