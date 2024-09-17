# Data generation pipeline for 'ground truth' simulations in PyBullet

using Dates
using Gen
using CSV
using DataFrames
using PyCall

bullet = pyimport("pybullet")
pybullet_data = pyimport("pybullet_data")

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

    #=
    frontWallID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.5, 0.01, 0.5])
    wall4 = bullet.createMultiBody(baseCollisionShapeIndex=frontWallID, basePosition=[0, -0.51, 0.51], baseOrientation=quaternion)
    bullet.changeDynamics(wall4, -1, mass=0, restitution=0.5)
    println("Front Wall", frontWallID)
    =#
    
    # Create and position ceiling
    ceilingID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[0.5, 0.5, 0.02])
    ceiling = bullet.createMultiBody(baseCollisionShapeIndex=ceilingID, basePosition=[0, 0, 1.01], baseOrientation=quaternion)
    bullet.changeDynamics(ceiling, -1, mass=0, restitution=0.5)
    println("Ceiling: ", ceilingID)

    # Set gravity
    bullet.setGravity(0, 0, -9.81)
    
end

# Stochastically sample the target object's initial kinematic state
function sample_target_state(elasticity)

    body = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.05, .05, .05])
    println("TargetID: ", body)

    position_range = [-0.4, 0.4]
    height_range = [0.1, 0.9]
    position = [uniform(position_range[1], position_range[2]), uniform(position_range[1], position_range[2]), uniform(height_range[1], height_range[2])]

    rotation_range = [0., 359.99]
    rotation = [uniform(rotation_range[1], rotation_range[2]), uniform(rotation_range[1], rotation_range[2]), uniform(rotation_range[1], rotation_range[2])]
    init_orientation = bullet.getQuaternionFromEuler(rotation)

    println(rotation)
    println(init_orientation)

    velocity_range = [-5., 5.]
    linear_velocity = [uniform(velocity_range[1], velocity_range[2]), uniform(velocity_range[1], velocity_range[2]), uniform(velocity_range[1], velocity_range[2])]

    target = bullet.createMultiBody(baseCollisionShapeIndex=body, basePosition=position, baseOrientation=init_orientation) 
    bullet.resetBaseVelocity(target, linearVelocity=linear_velocity)
    bullet.changeDynamics(target, -1, mass=1.0, restitution=elasticity)

    return target, rotation

end

function main()

    for i=3:29

        ela = div(i, 3)

        bullet.connect(bullet.GUI)::Int64
        bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
        bullet.resetDebugVisualizerCamera(3, 0, -3, [0, 0, 1])

        init_scene()

        elasticity = .1 * ela
        target, rotation = sample_target_state(elasticity)
        linear_velocity, angular_velocity = bullet.getBaseVelocity(target)
        init_position, orientation = bullet.getBasePositionAndOrientation(target)
        init_orientation = bullet.getEulerFromQuaternion(orientation)

        println(elasticity)
        println(orientation)
        println(init_orientation)

        state = DataFrame(
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
        
        # Create the initial state
        init_state = (
            Date = Date("2024-08-30"),
            SceneID = string("Cube_Ela", ela, "Var_", (i % 3 + 1)),
            Elasticity = elasticity,
            Variation = i % 3 + 1,
            PositionX = init_position[1],
            PositionY = init_position[2],
            PositionZ = init_position[3],
            RotationX = rotation[1],
            RotationY = rotation[2],
            RotationZ = rotation[3],
            VelocityX = linear_velocity[1],
            VelocityY = linear_velocity[2],
            VelocityZ = linear_velocity[3]
        )
        push!(state, init_state)
        println(state)

        # Truncate to 5 digits and write to CSV
        truncator(col, val) = trunc(val, digits=5)
        truncator(col, val::Int) = val
        CSV.write(string("Tests/BulletStimulus/Data/Cube_Ela", ela, "_Var", (i % 3 + 1), ".csv"), state)
        
        # Generate the trajectory
        trajectory_data = DataFrame(x=[], y=[], z=[])
        vel = linear_velocity
        while abs(vel[1]) > 0.1 || abs(vel[2]) > 0.1 || abs(vel[3]) > 0.1
            
            bullet.stepSimulation()
            
            position, orientation = bullet.getBasePositionAndOrientation(target)
            orientation = bullet.getEulerFromQuaternion(orientation)
            time_step = [position[1]; position[2]; position[3]]
            #println(time_step)
            push!(trajectory_data, time_step)

            aabbMin, aabbMax = bullet.getAABB(target)
            colliders = bullet.getOverlappingObjects(aabbMin, aabbMax)
            #println(length(colliders), ": ", colliders)

            vel, angular_velocity = bullet.getBaseVelocity(target)
            println(vel)

        end

        # Write to CSV
        CSV.write(string("Tests/BulletStimulus/Data/Cube_Ela", ela, "_Var", (i % 3 + 1), "_observed.csv"), trajectory_data, transform=truncator)

        bullet.disconnect()
        
    end

end

main()