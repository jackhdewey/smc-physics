# Generative Model
#
# Generates a rigid cube or sphere with specified initial position and velocity, within a larger cubic enclosure
# Simulates rigid body physics for specified number of time steps, recording a series of noisy position observations 
#
# DONE: Add noise to position
# TODO: Add noise to initial velocity
# CONSIDER: Add noise to orientation
# CONSIDER: Add additional noise to collisions
# CONSIDER: Change number of particles, number of rejuvenation moves
# CONSIDER: Changing prior for elasticity

using Gen
using PyCall
using PhySMC
using PhyBullet
using Accessors

bullet = pyimport("pybullet")
pybullet_data = pyimport("pybullet_data")

include("../Utilities/truncatednorm.jl")


# Sets initial scene configuration - in the future this could be inferred using an inverse graphics module
function init_scene()

    # Set gravity
    bullet.setGravity(0, 0, -9.81)

    # Create and position ground plane
    planeID = bullet.createCollisionShape(bullet.GEOM_PLANE)
    plane = bullet.createMultiBody(baseCollisionShapeIndex=planeID, basePosition=[0, 0, 0])
    bullet.changeDynamics(plane, -1, mass=0.0, restitution=0.5)

    # Create and position walls
    quaternion = bullet.getQuaternionFromEuler([0, 0, 0])

    wall1ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.01, 0.5, 0.5])
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
    
end

# Sets the initial state (shape, position, orientation, velocity, and optionally mass and restitution) of the target object
function init_target_state(sim::PhySim, shape::String, init_position::Vector{Float64}, init_velocity::Vector{Float64}, mass::Float64=1.0, restitution::Float64=0.9)
 
    # Select shape representation
    if (shape === "Cube")
        print("Initializing Cube\n")
        body = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.05, .05, .05])
    else
        print("Initializing Sphere\n")
        body = bullet.createCollisionShape(bullet.GEOM_SPHERE, radius=.1)
    end

    # Create, position, and orient target object
    init_orientation = bullet.getQuaternionFromEuler([0, 0, 0])
    target = bullet.createMultiBody(baseCollisionShapeIndex=body, basePosition=init_position, baseOrientation=init_orientation)    

    # Set initial velocity
    bullet.resetBaseVelocity(target, linearVelocity=init_velocity)

    # Set latent dynamic state
    bullet.changeDynamics(target, -1, mass=mass, restitution=restitution)
 
    # Store representation of cube's initial state
    init_state = BulletState(sim, [RigidBody(target)])

    return init_state
end

# Samples latent properties from their priors
@gen function sample_latents(l::RigidBodyLatents)

    res = {:restitution} ~ uniform(0, 1)
    
    # mass = {:mass} ~ gamma(1.2, 10.)

    return RigidBodyLatents(setproperties(l.data, restitution=res))
end

# Adds noise to kinematic state given transition uncertainty
# Current estimate as mean, variance either some constant or derived from average acceleration
@gen function resample_state(k::RigidBodyState)

    position = {:position} ~ broadcasted_normal(k.position, .01)

    #=
    orientation::Vector{3, Float64} = bullet.getEulerFromQuaternion(k.orientation)
    orientation = {:orientation} ~ broadcasted_normal(orientation, 0.1)
    orientation = bullet.getQuaternionFromEuler(orientation)
    =#

    return setproperties(k, position=position)
end

# Adds measurement noise to estimated position
@gen function generate_observation(k::RigidBodyState)

    obs = {:position} ~ broadcasted_normal(k.position, 0.2)

    return obs
end

# Given an input state, samples an observation and generates the next state
# Before calling PhySMC.step, perturb position / velocity / orientation and store in new state
@gen function kernel(t::Int, current_state::BulletState, sim::BulletSim)

    # Applies system noise to position, orientation, and velocity
    noisy_kinematics = {:state} ~ Gen.Map(resample_state)(current_state.kinematics)
    current_state = setproperties(current_state, kinematics = noisy_kinematics)

    # Applies observation noise to x, y, and z position
    {:observation} ~ Gen.Map(generate_observation)(current_state.kinematics)

    # Synchronizes state, then use Bullet to generate next state
    next_state::BulletState = PhySMC.step(sim, current_state)

    return next_state
end

#= 
TODO: Perturb initial kinematic state before calling Gen.Unfold
    - Ground truth as mean, variance derived from empirical distribution of data 
=#
# Given an initial state, samples latents from their priors then runs a complete forward simulation
@gen function generate_trajectory(sim::BulletSim, init_state::BulletState, T::Int)

    # Sample the target object's restitution
    latents = {:latents} ~ Gen.Map(sample_latents)(init_state.latents)
    init_state = setproperties(init_state; latents=latents)

    # Simulate T time steps
    states = {:trajectory} ~ Gen.Unfold(kernel)(T, init_state, sim)

    return states
end

# Generates a constrained trace with the specified arguments and elasticity
function generate_constrained(args::Tuple, restitution::Float64)
  
    gt_constraints = choicemap((:latents => 1 => :restitution, restitution))
    ground_truth = first(generate(generate_trajectory, args, gt_constraints))
        
    return ground_truth
end
