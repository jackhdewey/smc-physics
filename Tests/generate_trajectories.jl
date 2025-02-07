# Data generation pipeline for 'ground truth' simulations in PyBullet

using Dates
using Gen
using PyCall
using PhySMC
using PhyBullet
using CSV
using DataFrames

bullet = pyimport("pybullet")
pybullet_data = pyimport("pybullet_data")

obj_type::String = "Cube"
debug_viz::Bool = false

# Sets initial scene configuration (in the future this could be inferred using an inverse graphics module)
function init_scene()

    # Create and position ground plane
    planeID = bullet.createCollisionShape(bullet.GEOM_PLANE)
    plane = bullet.createMultiBody(baseCollisionShapeIndex=planeID, basePosition=[0, 0, 0])
    bullet.changeDynamics(plane, -1, mass=0, restitution=0.5)
    println("Ground Plane: ", planeID)

    # Create and position walls
    quaternion = bullet.getQuaternionFromEuler([0, 0, 0])

    leftWallID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.01, 0.5, 0.5])
    wall1 = bullet.createMultiBody(baseCollisionShapeIndex=leftWallID, basePosition=[-0.51, 0, 0.51], baseOrientation=quaternion)
    bullet.changeDynamics(wall1, -1, mass=0, restitution=0.5)
    println("Left Wall: ", leftWallID)

    rightWallID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.01, 0.5, 0.5])
    wall3 = bullet.createMultiBody(baseCollisionShapeIndex=rightWallID, basePosition=[0.51, 0, 0.51], baseOrientation=quaternion)
    bullet.changeDynamics(wall3, -1, mass=0, restitution=0.5)
    println("Right Wall: ", rightWallID)

    backWallID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.5, 0.01, 0.5])
    wall2 = bullet.createMultiBody(baseCollisionShapeIndex=backWallID, basePosition=[0, 0.51, 0.51], baseOrientation=quaternion)
    bullet.changeDynamics(wall2, -1, mass=0, restitution=0.5)
    println("Back Wall: ", backWallID)

    if !debug_viz
        frontWallID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.5, 0.01, 0.5])
        wall4 = bullet.createMultiBody(baseCollisionShapeIndex=frontWallID, basePosition=[0, -0.51, 0.51], baseOrientation=quaternion)
        bullet.changeDynamics(wall4, -1, mass=0, restitution=0.5)
        println("Front Wall", frontWallID)
    end
    
    # Create and position ceiling
    ceilingID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.5, 0.5, 0.02])
    ceiling = bullet.createMultiBody(baseCollisionShapeIndex=ceilingID, basePosition=[0, 0, 1.01], baseOrientation=quaternion)
    bullet.changeDynamics(ceiling, -1, mass=0, restitution=0.5)
    println("Ceiling: ", ceilingID)

    # Set gravity
    bullet.setGravity(0, 0, -9.81)
    
end

# Stochastically sample the target object's initial kinematic state
function sample_target_state(elasticity, sim)

    if obj_type == "Cube"
        body = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.05, .05, .05])
    else 
        body = bullet.createCollisionShape(bullet.GEOM_SPHERE, radius=.05)
    end
    println("TargetID: ", body)

    position_range = [-0.4, 0.4]
    height_range = [0.1, 0.9]
    position = [uniform(position_range[1], position_range[2]), 
                uniform(position_range[1], position_range[2]), 
                uniform(height_range[1], height_range[2])]

    rotation_range = [0., 2pi]
    orientation_euler = [uniform(rotation_range[1], rotation_range[2]), 
                uniform(rotation_range[1], rotation_range[2]), 
                uniform(rotation_range[1], rotation_range[2])]
    orientation_quat = bullet.getQuaternionFromEuler(orientation_euler)

    println("Orientation: ", orientation_euler)
    println("Orientation: ", orientation_quat)

    velocity_range = [-5., 5.]
    linear_velocity = [ uniform(velocity_range[1], velocity_range[2]), 
                        uniform(velocity_range[1], velocity_range[2]), 
                        uniform(velocity_range[1], velocity_range[2])]

    target = bullet.createMultiBody(baseCollisionShapeIndex=body, basePosition=position, baseOrientation=orientation_quat) 
    bullet.resetBaseVelocity(target, linearVelocity=linear_velocity)
    bullet.changeDynamics(target, -1, mass=1.0, restitution=elasticity)

    init_state = BulletState(sim, [RigidBody(target)])

    return target, init_state

end

function main()

    for i=3:29

        # Establish bullet server
        if debug_viz
            client = bullet.connect(bullet.GUI)::Int64
        else
            client = bullet.connect(bullet.DIRECT)::Int64
        end
        bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
        bullet.resetDebugVisualizerCamera(3, 0, -3, [0, 0, 1])

        # Initilize PhyBullet
        sim = BulletSim(step_dur=1/60; client=client)

        # Set up container room
        init_scene()

        # Set elasticity, initial state
        ela = div(i, 3)
        elasticity = .1 * ela
        target, init_state = sample_target_state(elasticity, sim)

        # Print initial state
        linear_velocity, angular_velocity = bullet.getBaseVelocity(target)
        println("Elasticity: ", elasticity)
        init_position, orientation = bullet.getBasePositionAndOrientation(target)
        println("Orientation (Quaternion): ", orientation)
        init_orientation = bullet.getEulerFromQuaternion(orientation)
        println("Orientation (Euler): ", init_orientation)

        # Setup and initialize state representation
        state_data = DataFrame(
            Date=Date[], 
            SceneID=String[], 
            Elasticity=Float64[], 
            Variation=Int[], 
            PositionX=Float64[], 
            PositionY=Float64[], 
            PositionZ=Float64[], 
            RotationX=Float64[], 
            RotationY=Float64[], 
            RotationZ=Float64[], 
            VelocityX=Float64[], 
            VelocityY=Float64[], 
            VelocityZ=Float64[]
        )
        init_state_data = (
            Date = Date("2024-08-30"),
            SceneID = string("Cube_Ela", ela, "Var_", (i % 3 + 1)),
            Elasticity = elasticity,
            Variation = i % 3 + 1,
            PositionX = init_position[1],
            PositionY = init_position[2],
            PositionZ = init_position[3],
            RotationX = init_orientation[1],
            RotationY = init_orientation[2],
            RotationZ = init_orientation[3],
            VelocityX = linear_velocity[1],
            VelocityY = linear_velocity[2],
            VelocityZ = linear_velocity[3]
        )
        push!(state_data, init_state_data)

        # Truncate to 5 digits and write to CSV
        truncator(col, val) = trunc(val, digits=5)
        truncator(col, val::Int) = val
        CSV.write(string("Data/BulletStimulus/", obj_type, "/", obj_type, "_Ela", ela, "_Var", (i % 3 + 1), ".csv"), state_data)
        
        # Generate the trajectory
        t = 0
        trajectory_data = DataFrame(x=[], y=[], z=[])
        vel = linear_velocity
        state::BulletState = init_state
        while t < 120 && (abs(vel[1]) > 0.1 || abs(vel[2]) > 0.1 || abs(vel[3]) > 0.1)
            
            state = PhySMC.step(sim, state)
            #bullet.stepSimulation()
     
            position, orientation = bullet.getBasePositionAndOrientation(target)
            println("Position: ", position)

            orientation = bullet.getEulerFromQuaternion(orientation)
            println("Orientation: ", orientation)

            push!(trajectory_data, position)

            aabbMin, aabbMax = bullet.getAABB(target)
            colliders = bullet.getOverlappingObjects(aabbMin, aabbMax)
            #println(length(colliders), ": ", colliders)

            vel, angular_velocity = bullet.getBaseVelocity(target)
            #println(vel)

            t += 1
        end

        # Write to CSV
        CSV.write(string("Data/BulletStimulus/New/", obj_type, "/", obj_type, "_Ela", ela, "_Var", (i % 3 + 1), "_observed.csv"), trajectory_data, transform=truncator)
        
        bullet.disconnect()
        
    end

end

main()