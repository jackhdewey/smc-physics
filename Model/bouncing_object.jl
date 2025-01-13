# Generative Model
#
# Generates a PyBullet enviorment featuring a rigid cube or sphere (with variable initial position and velocity) within a larger cubic enclosure
# Simulates rigid body physics for a specified number of time steps, returning a series of (noisy) position observations 
#
# CONSIDER: Add noise to initial orientation
# CONSIDER: Add additional noise at collisions

using Gen
using PyCall
using PhySMC
using PhyBullet
using Accessors

bullet = pyimport("pybullet")
pybullet_data = pyimport("pybullet_data")

include("../args.jl")
include("../Utilities/truncatednorm.jl")


model_args = Args()

# Sets initial scene configuration (in the future this could be inferred using an inverse graphics module)
function init_scene(debug_viz::Bool)

    # Create and position ground plane
    planeID = bullet.createCollisionShape(bullet.GEOM_PLANE)
    plane = bullet.createMultiBody(baseCollisionShapeIndex=planeID, basePosition=[0, 0, 0])
    bullet.changeDynamics(plane, -1, mass=0.0, restitution=0.5)

    # Create and position walls
    quaternion = bullet.getQuaternionFromEuler([0, 0, 0])

    wall1ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.01, 0.5, 0.5])
    wall1 = bullet.createMultiBody(baseCollisionShapeIndex=wall1ID, basePosition=[-0.51, 0, 0.51], baseOrientation=quaternion)
    bullet.changeDynamics(wall1, -1, mass=0.0, restitution=0.5)

    wall2ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.5, 0.01, 0.5])
    wall2 = bullet.createMultiBody(baseCollisionShapeIndex=wall2ID, basePosition=[0, 0.51, 0.51], baseOrientation=quaternion)
    bullet.changeDynamics(wall2, -1, mass=0.0, restitution=0.5)

    wall3ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.01, 0.5, 0.5])
    wall3 = bullet.createMultiBody(baseCollisionShapeIndex=wall3ID, basePosition=[0.51, 0, 0.51], baseOrientation=quaternion)
    bullet.changeDynamics(wall3, -1, mass=0.0, restitution=0.5)

    if !debug_viz
        wall4ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.5, 0.01, 0.5])
        wall4 = bullet.createMultiBody(baseCollisionShapeIndex=wall4ID, basePosition=[0, -0.51, 0.51], baseOrientation=quaternion)
        bullet.changeDynamics(wall4, -1, mass=0.0, restitution=0.5)
    end
    
    # Create and position ceiling
    ceilingID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.5, 0.5, 0.01])
    ceiling = bullet.createMultiBody(baseCollisionShapeIndex=ceilingID, basePosition=[0, 0, 1.01], baseOrientation=quaternion)
    bullet.changeDynamics(ceiling, -1, mass=0.0, restitution=0.5)

    # Set gravity
    bullet.setGravity(0, 0, -9.81)
    
end

# Sets target object's initial state (shape, position, orientation, velocity, and optionally mass and restitution)
function init_target_state(sim::PhySim, shape::String, init_position::Vector{Float64}, init_velocity::Vector{Float64}, mass::Float64=1.0, restitution::Float64=0.9)
 
    # Select shape representation
    if (shape === "Cube")
        print("Initializing Cube\n")
        body = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.05, .05, .05])
    else
        print("Initializing Sphere\n")
        body = bullet.createCollisionShape(bullet.GEOM_SPHERE, radius=.05)
    end

    # Create targets object with initial position and orientation
    init_orientation = bullet.getQuaternionFromEuler([0, 0, 0])
    target = bullet.createMultiBody(baseCollisionShapeIndex=body, basePosition=init_position, baseOrientation=init_orientation)    

    # Set initial velocity and latent dynamic properties
    bullet.resetBaseVelocity(target, linearVelocity=init_velocity)
    bullet.changeDynamics(target, -1, mass=mass, restitution=restitution)
 
    # Create PhyBullet representation of target's initial state
    init_state = BulletState(sim, [RigidBody(target)])

    return init_state
end

# Samples the target object's latent dynamic properties from their priors
@gen function sample_latents(l::RigidBodyLatents)

    res = {:restitution} ~ uniform(0, 1)

    # mass = {:mass} ~ gamma(1.2, 10.)

    return RigidBodyLatents(setproperties(l.data, restitution=res))
end

# Samples the target object's initial velocity
@gen function sample_init_state(k::RigidBodyState)

    velocity = {:init_velocity} ~ broadcasted_normal(k.linear_vel, model_args.init_vel_noise)

    return setproperties(k, linear_vel=velocity)
end

# Adds noise to position to capture uncertainty in transition dynamics
#   - Current estimate as mean, variance constant 
#   - Alternatives: ground truth as mean; derive variance from average acceleration (empirical distribution of data)
@gen function resample_state(k::RigidBodyState)

    position = {:position} ~ broadcasted_normal(k.position, model_args.transition_noise)
    
    #=
    orientation::Vector{3, Float64} = bullet.getEulerFromQuaternion(k.orientation)
    orientation = {:orientation} ~ broadcasted_normal(orientation, 0.1)
    orientation = bullet.getQuaternionFromEuler(orientation)
    =#

    return setproperties(k, position=position)
end

# Adds measurement noise to position to capture uncertainty in observation
@gen function sample_observation(k::RigidBodyState)

    obs = {:position} ~ broadcasted_normal(k.position, model_args.observation_noise)

    return obs
end

# Core transition and observation kernel - perturbs the current state, samples an observation, and generates the next state using simulator
# Before calling PhySMC.step, perturb position / velocity / orientation and store in new state
@gen function kernel(t::Int, current_state::BulletState, sim::BulletSim)

    # Applies system noise to predicted position
    noisy_kinematics = {:state} ~ Gen.Map(resample_state)(current_state.kinematics)
    current_state = setproperties(current_state, kinematics=noisy_kinematics)

    # Applies observation noise to predicted position generate predicted observation
    {:observation} ~ Gen.Map(sample_observation)(current_state.kinematics)

    # Synchronizes state then calls Bullet to generate next state
    next_state::BulletState = PhySMC.step(sim, current_state)

    return next_state
end

# Given an initial state, samples latents from their priors then runs a stochastic forward simulation
@gen function generate_trajectory(sim::BulletSim, init_state::BulletState, T::Int)

    # Sample the object's restitution (elasticity)
    latents = {:latents} ~ Gen.Map(sample_latents)(init_state.latents)
    init_state = setproperties(init_state; latents=latents)

    # Sample the object's initial velocity
    #kinematics = { :init_state } ~ Gen.Map(sample_init_state)(init_state.kinematics)
    #init_state = setproperties(init_state; kinematics=kinematics)

    # Simulate T time steps
    states = {:trajectory} ~ Gen.Unfold(kernel)(T, init_state, sim)

    return states
end

# Generates a trace with the specified elasticity
function generate_constrained(sim_args::Tuple, restitution::Float64)
  
    gt_constraints = choicemap((:latents => 1 => :restitution, restitution))
    ground_truth = first(generate(generate_trajectory, sim_args, gt_constraints))
        
    return ground_truth
end
