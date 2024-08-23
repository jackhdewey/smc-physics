# Data generation pipeline for 'ground truth' simulations in PyBullet

using CSV
using DataFrames

include("../Model/bouncing_object.jl")

function init_target_state(elasticity)

    position_range = [-0.45, 0.45]
    height_range = [0.05, 0.95]
    rotation_range = [0., 359.99]
    velocity_range = [-2., 2.]

    position = [uniform(position_range[1], position_range[2]), uniform(position_range[1], position_range[2]), uniform(height_range[1], height_range[2])]
    rotation = [uniform(rotation_range[1], rotation_range[2]), uniform(rotation_range[1], rotation_range[2]), uniform(rotation_range[1], rotation_range[2])]
    velocity = [uniform(velocity_range[1], velocity_range[2]), uniform(velocity_range[1], velocity_range[2]), uniform(velocity_range[1], velocity_range[2])]

    body = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.05, .05, .05])
    init_orientation = bullet.getQuaternionFromEuler(rotation)
    target = bullet.createMultiBody(baseCollisionShapeIndex=body, basePosition=position, baseOrientation=init_orientation) 
    bullet.resetBaseVelocity(target, linearVelocity=velocity)
    bullet.changeDynamics(target, -1, restitution=elasticity)

    return target

end

function main()

    #sim = BulletSim(step_dur=1/30; client=client)
    #for i=3:3

    client = bullet.connect(bullet.GUI)::Int64
    bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
    bullet.resetDebugVisualizerCamera(3, 0, -3, [0, 0, 2])
    init_scene()

    elasticity = 1
    println(elasticity)
    target = init_target_state(elasticity)
    velocity = bullet.getBaseVelocity(target)
    vel = velocity[1]
    trajectory_data = DataFrame(elasticity=[], x=[], y=[], z=[], x0=[], y0=[], z0=[])

    while vel[1] > 0. || vel[2] > 0. || vel[3] > 0.
        bullet.stepSimulation()
        position, orientation = bullet.getBasePositionAndOrientation(target)
        orientation = bullet.getEulerFromQuaternion(orientation)
        #pos = convert(Vector, position)
        #rot = convert(Vector, orientation)
        time_step = [elasticity; position[1]; position[2]; position[3]; orientation[1]; orientation[2]; orientation[3]]
        println(time_step)
        push!(trajectory_data, time_step)
        velocity = bullet.getBaseVelocity(target)
        vel = velocity[1]
    end

    # Truncate to 5 digits
    truncator(col, val) = trunc(val, digits=5)
    truncator(col, val::Int) = val

    CSV.write(string("test", i, ".csv"), trajectory_data, transform=truncator)

    bullet.disconnect()
        
    #end

end

main()