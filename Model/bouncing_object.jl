# Generative Model
#
# Generates a rigid cube or sphere with specified initial position and velocity, within a larger cubic enclosure
# Simulates rigid body physics for specified number of time steps, recording a series of noisy position observations 
#
# TODO: Add noise to initial position, orientation, and velocity
# TODO: Add noise to velocity at each time step
# TODO: Add additional noise to collisions
# TODO: Change number of particles,  number of reuvenation moves
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

# Samples an initial estimate of latent properties
@gen function sample_latents(latents::RigidBodyLatents)

    # mass = {:mass} ~ gamma(1.2, 10.)
    res = {:restitution} ~ uniform(0, 1)
    new_latents = RigidBodyLatents(setproperties(latents.data, restitution=res))

    return new_latents
end

# Sets the initial state (shape, position, orientation, velocity, and optionally mass and restitution) of the target object
@gen function init_target_state(sim::PhySim, shape::String, init_position::Vector{Float64}, init_velocity::Vector{Float64}, mass::Float64=1.0, restitution::Float64=0.9)
 
    # Select shape representation
    if (shape === "Cube")
        print("Initializing Cube\n")
        body = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.05, .05, .05])
    else
        print("Initializing Sphere\n")
        body = bullet.createCollisionShape(bullet.GEOM_SPHERE, radius=.1)
    end

    # Sample initial kinematic state - position, orientation, and velocity

    # init_velocity = {:init_velocity} ~
    # startPosition = {:init_position} ~
    # startOrientation = {:init_orientation} ~
    
    startOrientation = bullet.getQuaternionFromEuler([0, 0, 1])
 
    # Create, position, and orient target object
    target = bullet.createMultiBody(baseCollisionShapeIndex=body, basePosition=init_position, baseOrientation=startOrientation)    

    bullet.resetBaseVelocity(target, linearVelocity=init_velocity)

    # Set latent dynamic state
    bullet.changeDynamics(target, -1, mass=mass, restitution=restitution)
 
    # Store representation of cube's initial state
    init_state = BulletState(sim, [RigidBody(target)])

    return init_state
end

# Adds measurement noise to ground truth position
# TODO: Increase variance
@gen function generate_observation(k::RigidBodyState)

    obs = {:position} ~ broadcasted_normal(k.position, 0.1)

    return obs
end

#=
Option:
    - Before calling PhySMC.step, perturb position / velocity / orientation and store in new state
        - Current estimate as mean, variance either some constant or derived from average acceleration
=#

# Given an input state, samples an observation and generates the next state
@gen function kernel(t::Int, current_state::BulletState, sim::BulletSim)

    # Applies noise to x, y, and z positions
    {:observation} ~ Gen.Map(generate_observation)(current_state.kinematics)

    # Synchronizes state, then use Bullet to generate next state
    next_state::BulletState = PhySMC.step(sim, current_state)

    return next_state
end

#= 
Option:
    - Perturb initial kinematic state before calling Gen.Unfold
        - Ground truth as mean, variance derived from empirical distribution of data 
=#

# Given an initial state, samples latents from their priors then runs a complete forward simulation
@gen function generate_trajectory(sim::BulletSim, init_state::BulletState, T::Int)

    # Sample the target object's latents - i.e. restitution
    latents = {:latents} ~ Gen.Map(sample_latents)(init_state.latents)

    # Update initial state with sampled latent values
    init_state = Accessors.setproperties(init_state; latents=latents)

    # Simulate T time steps
    states = {:trajectory} ~ Gen.Unfold(kernel)(T, init_state, sim)

    return states
end

# Generates a constrained trace with the specified arguments and elasticity
function generate_ground_truth(args::Tuple, restitution::Float64)
  
    gt_constraints = choicemap((:latents => 1 => :restitution, restitution))
    ground_truth = first(generate(generate_trajectory, args, gt_constraints))
        
    return ground_truth
end
