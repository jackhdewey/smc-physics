# Functions for testing the behavior of PyBullet solver under different initial conditions

# TODO: Compare RealFlow and PyBullet bounce heights for different elasticities in very simple example (dropping a sphere)
#   * Drop sphere from 1m, vary elasticity from 0.0-1.0 with 0.1 increments, save to .csv 

# TODO: Compare RealFlow and PyBullet bounce directions for different initial orientations
#   * Rotate a cube to ~10 angles, then do the same test

import pybullet as p
import time
import numpy
import csv

# Create a scene with a ground plane and four walls enclosing a cube / sphere dropped from 1m 
def init_scene():

    # Create and position ground plane
    planeID = p.createCollisionShape(p.GEOM_PLANE)
    plane = p.createMultiBody(baseCollisionShapeIndex=planeID, basePosition=[0, 0, 0])
    p.changeDynamics(plane, -1, restitution=0.3)

    # Create and position walls
    wall1ID = p.createCollisionShape(p.GEOM_BOX, halfExtents=[0.01, .5, .5])
    quaternion = p.getQuaternionFromEuler([0, 0, 0])
    wall1 = p.createMultiBody(baseCollisionShapeIndex=wall1ID, basePosition=[.5, 0, .5], baseOrientation=quaternion)
    p.changeDynamics(wall1, -1, restitution=0.3)

    wall2ID = p.createCollisionShape(p.GEOM_BOX, halfExtents=[0.01, .5, .5])
    wall2 = p.createMultiBody(baseCollisionShapeIndex=wall2ID, basePosition=[-.5, 0, .5], baseOrientation=quaternion)
    p.changeDynamics(wall2, -1, restitution=0.3)

    wall3ID = p.createCollisionShape(p.GEOM_BOX, halfExtents=[.5, 0.01, .5])
    wall3 = p.createMultiBody(baseCollisionShapeIndex=wall3ID, basePosition=[0, .5, .5], baseOrientation=quaternion)
    p.changeDynamics(wall3, -1, restitution=0.3)

    '''
    wall4ID = p.createCollisionShape(p.GEOM_BOX, halfExtents=[.5, 0.01, .5])
    wall4 = p.createMultiBody(baseCollisionShapeIndex=wall4ID, basePosition=[0, -.5, .5], baseOrientation=quaternion)
    p.changeDynamics(wall4, -1, restitution=0.3)
    '''
    
    rigidBodies = []
    init_locations = [[-0.6, 0, 1], [-0.2, 0, 1], [0.2, 0, 1], [0.6, 0, 1]]

    # Create cube collision shape
    cubeShape = p.createCollisionShape(p.GEOM_BOX, halfExtents=[.05, .05, .05])

    """ # Position and orient cube
    startPosCube = [0, 0, 1]
    startOrientationCube = p.getQuaternionFromEuler([0, 0, 0])
    cube = p.createMultiBody(baseCollisionShapeIndex=cubeShape, basePosition=startPosCube, baseOrientation=startOrientationCube)

    # Set dynamics 
    p.setGravity(0, 0, -9.81)
    p.changeDynamics(cube, -1, mass=1.0, restitution=0.9)

    rigidBodies.append(cube) """


    # Create sphere collision shape
    sphereShape = p.createCollisionShape(p.GEOM_SPHERE, radius=.062)

    # Position and orient sphere
    startPosSphere = [0, -1, 1]
    startOrientationSphere = p.getQuaternionFromEuler([0, 0, 0])
    sphere = p.createMultiBody(baseCollisionShapeIndex=sphereShape, basePosition=startPosSphere, baseOrientation=startOrientationSphere)

    p.changeDynamics(sphere, -1, mass=1.0, restitution=0.9)

    rigidBodies.append(sphere)

    '''
    # Position and orient second cube
    startPosCube = cube_locations[1]
    cube1 = p.createMultiBody(baseCollisionShapeIndex=cubeShape, basePosition=startPosCube,
                              baseOrientation=startOrientationCube)
    p.changeDynamics(cube1, -1, mass=1.0, restitution=0.7)
    cubeBodies.append(cube1)
    
    # Create and position soft body cubes
    cube1 = p.loadSoftBody("soft_cube__sf.obj", simFileName="soft_cube.vtk", basePosition=cube_locations[0],
                            mass=4, useNeoHookean=1, NeoHookeanMu=60, NeoHookeanLambda=600, NeoHookeanDamping=0.001,
                            collisionMargin=0.006, useSelfCollision=1, frictionCoeff=0.5, repulsionStiffness=800)
    cubeBodies.append(cube1)

    cube2 = p.loadSoftBody("soft_cube__sf.obj", simFileName="soft_cube.vtk", basePosition=cube_locations[1],
                            mass=4, useNeoHookean=1, NeoHookeanMu=280, NeoHookeanLambda=600, NeoHookeanDamping=0.001,
                            collisionMargin=0.006, useSelfCollision=1, frictionCoeff=0.5, repulsionStiffness=700)
    cubeBodies.append(cube2)
    '''

    return rigidBodies

############### TESTS ################

# Run a basic simulation
def run_sample(cubeBodies, cube_positions):

    for _ in range(10000):
        p.stepSimulation()
        current_position = p.getBasePositionAndOrientation(cubeBodies[0])[0]
        cube_positions.append(current_position)
        time.sleep(1. / 3640.)

    return cube_positions

# Run ten simulations with elasticity progressively increased from 0.0 to 1.0
def test_elasticity(object, shape):

    # Open a .csv file
    file = "Tests/test_elasticity_" + shape + ".csv"
    csvfile = open(file, 'w')
    writer = csv.writer(csvfile)

    # Write the fields to the first row
    fields = ['Elasticity', 'X', 'Y', 'Z']
    writer.writerow(fields)

    for i in range(0, 10):

        # Procedurally vary elasticity
        elasticity = 0.1 * i
        p.changeDynamics(object, -1, restitution=elasticity)

        # Concatenate initial position and velocity
        initial_position = [0.0, 0.0, 1.0]
        initial_position.insert(0, elasticity)
        writer.writerow(initial_position)

        # Simulate 200 time steps and record trajectory
        for _ in range(0, 200):

            p.stepSimulation()
            
            position = p.getBasePositionAndOrientation(object)[0]
            position = list(position)
            position.insert(0, elasticity)
            writer.writerow(position)       

            time.sleep(1. / 3640.) 

        # Reset to starting position / orientation
        p.resetBasePositionAndOrientation(object, [0, 0, 1], p.getQuaternionFromEuler([0, 0, 0]))

    return

# Run ten simulations with orientation progressively increased from 0.0 to pi/2
def test_orientation(object, shape):

    file = file = "Tests/test_orientation_" + shape + ".csv"
    csvfile = open(file, 'w')
    writer = csv.writer(csvfile)

    fields = ['Theta Y', 'X', 'Y', 'Z', 'Theta X', 'Theta Y', 'Theta Z']
    writer.writerow(fields)

    # Divide range from 0 to pi / 2 into ten intervals
    intervals = numpy.linspace(0, numpy.pi / 2, 10)

    # Run ten simulations 
    for i in range(0, 10):

        initial_position = [0.0, 0.0, 1.0]
        theta_y_init = round(intervals[i], 3)
        initial_orienation = [0.0, theta_y_init, 0.0]

        # Record initial position and orientation
        initial_position.insert(0, theta_y_init)
        initial_position.extend(initial_orienation)
        writer.writerow(initial_position)

        # Set to correct position and orientation in simulation
        p.resetBasePositionAndOrientation(object, [0, 0, 1], p.getQuaternionFromEuler([0, theta_y_init, 0]))

        # Simulate 200 time steps and record trajectory
        for _ in range(0, 200):

            p.stepSimulation()

            # Extract current position and orientation
            position, quaternion = p.getBasePositionAndOrientation(object)
            orientation = p.getEulerFromQuaternion(quaternion)
            position = list(position)
            position.insert(0, theta_y_init)
            position.extend(orientation)
            
            # Round values to three decimals points
            for i in range(len(position)):
                position[i] = round(position[i], 3)

            # Write to .csv
            writer.writerow(position)        

            time.sleep(1. / 3640.)

    return
