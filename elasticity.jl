# Generative Model
#
# Generates an elastic cube or sphere with specified position and velocity within a larger cubic enclosure
# Simulates rigid body physics for specified number of time steps, recording a series of noisy position observations
# Implemented as a deterministic simulation over uncertain properties
#
# CONSIDER: Changing prior for elasticity
# CONSIDER: Adding noise over velocity
# CONSIDER: Handling collisions separately

using Gen
using PyCall
using PhySMC
using PhyBullet
using Accessors

bullet = pyimport("pybullet")
pybullet_data = pyimport("pybullet_data")
include("Utilities/truncatednorm.jl")


# Sets initial scene configuration
# In the future this will be an inverse graphics module
function init_scene()

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

    bullet.setGravity(0, 0, -10)

end

# Sets the initial state (shape, position, velocity, and optionally mass and restitution of the target object)
function init_target_state(sim::PhySim, shape::String, init_position::Vector{Float64}, init_velocity::Vector{Float64}, mass::Float64=1.0, restitution::Float64=0.9)

     # Set target object initial position and orientation
     startPosition = init_position
     startOrientation = bullet.getQuaternionFromEuler([0, 0, 1])
 
     # Select shape representation
     if (shape === "Cube")
        print("Initializing Cube\n")
        body = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.05, .05, .05])
     else
        print("Initializing Sphere\n")
        body = bullet.createCollisionShape(bullet.GEOM_SPHERE, radius=.1)
     end
 
     # Create, position, and orient target object
     target = bullet.createMultiBody(baseCollisionShapeIndex=body, basePosition=startPosition, baseOrientation=startOrientation)
 
     # Set dynamic properties and initial kinematic state
     bullet.changeDynamics(target, -1, mass=mass, restitution=restitution)
     bullet.resetBaseVelocity(target, linearVelocity=init_velocity)
 
     # Store representation of cube's initial state
     init_state = BulletState(sim, [RigidBody(target)])

     return init_state
end

# Samples an initial estimate of latent properties
@gen function sample_latents(latents::RigidBodyLatents)
    # mass = {:mass} ~ gamma(1.2, 10.)
    res = {:restitution} ~ uniform(0, 1)
    new_latents = RigidBodyLatents(setproperties(latents.data, restitution=res))
    return new_latents
end

# Adds measurement noise to ground truth position
@gen function generate_observation(k::RigidBodyState)
    obs = {:position} ~ broadcasted_normal(k.position, 0.1)
    return obs
end

# Given an input state, samples an observation, syncs Bullet simulator, and generates the next state
@gen function kernel(t::Int, current_state::BulletState, sim::BulletSim)

    # Applies noise to x, y, and z positions
    {:observation} ~ Gen.Map(generate_observation)(current_state.kinematics)

    # Synchronizes state, then use Bullet to generate next state
    next_state::BulletState = PhySMC.step(sim, current_state)

    return next_state
end

# Given an initial state, samples latents from their priors then runs a complete forward simulation
@gen function generate_trajectory(sim::BulletSim, init_state::BulletState, T::Int)

    # Sample values for the target object's latents - i.e. a restitution
    latents = {:latents} ~ Gen.Map(sample_latents)(init_state.latents)

    # Update initial state with sampled latent values
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
