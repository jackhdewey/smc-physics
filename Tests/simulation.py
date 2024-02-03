# Functions for testing the behavior of PyBullet solver under different initial conditions

# DONE: Set elasticity of floor to 0.3

# DONE: Compare RealFlow and PyBullet bounce heights for different elasticities in very simple example (dropping a sphere)
#   * Drop sphere from 1m, vary elasticity from 0.0-1.0 with 0.1 increments, save to .csv 

# DONE: Compare RealFlow and PyBullet bounce directions for different initial orientations
#   * Rotate a cube to ~10 angles, then do the same test

import pybullet as p
import time
import numpy
import csv

def init_scene():

    p.setGravity(0, 0, -10)

    # Create and position ground plane
    planeID = p.createCollisionShape(p.GEOM_PLANE)
    plane = p.createMultiBody(baseCollisionShapeIndex=planeID, basePosition=[0, 0, 0])
    p.changeDynamics(plane, -1, restitution=0.3)

    # Create and position walls
    plane2ID = p.createCollisionShape(p.GEOM_BOX, halfExtents=[0.01, 1, 1])
    quaternion = p.getQuaternionFromEuler([0, 0, 0])
    plane2 = p.createMultiBody(baseCollisionShapeIndex=plane2ID, basePosition=[1, 0, 1], baseOrientation=quaternion)
    p.changeDynamics(plane2, -1, restitution=0.3)

    plane3ID = p.createCollisionShape(p.GEOM_BOX, halfExtents=[0.01, 1, 1])
    plane3 = p.createMultiBody(baseCollisionShapeIndex=plane3ID, basePosition=[-1, 0, 1], baseOrientation=quaternion)
    p.changeDynamics(plane3, -1, restitution=0.3)

    plane4ID = p.createCollisionShape(p.GEOM_BOX, halfExtents=[1, 0.01, 1])
    plane4 = p.createMultiBody(baseCollisionShapeIndex=plane4ID, basePosition=[0, 1, 1], baseOrientation=quaternion)
    p.changeDynamics(plane4, -1, restitution=0.3)

    # Create base cube shape
    cubeShape = p.createCollisionShape(p.GEOM_BOX, halfExtents=[.1, .1, .1])
    cube_locations = [[-0.6, 0, 3], [-0.2, 0, 3], [0.2, 0, 3], [0.6, 0, 3]]
    startOrientationCube = p.getQuaternionFromEuler([0, 0, 0])
    cubeBodies = []
    
    # Position and orient first cube
    startPosCube = [0, 0, 1]
    cube = p.createMultiBody(baseCollisionShapeIndex=cubeShape, basePosition=startPosCube,
                             baseOrientation=startOrientationCube)
    p.changeDynamics(cube, -1, mass=1.0, restitution=0.9)
    cubeBodies.append(cube)

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

    return cubeBodies, cube_locations


def run_sample(cubeBodies, cube_positions):

    for i in range(10000):
        p.stepSimulation()
        current_position = p.getBasePositionAndOrientation(cubeBodies[0])[0]
        cube_positions.append(current_position)
        time.sleep(1. / 3640.)

    return


def test_elasticity(cube):

    fields = ['Elasticity', 'X', 'Y', 'Z']

    file = "Tests/test_elasticity.csv"
    csvfile = open(file, 'w')
    writer = csv.writer(csvfile)
    writer.writerow(fields)

    for i in range(0, 10):

        elasticity = 0.1 * i

        initial_position = [0.0, 0.0, 1.0]
        initial_position.insert(0, elasticity)
        writer.writerow(initial_position)

        for t in range(0, 200):

            p.stepSimulation()
            
            position = p.getBasePositionAndOrientation(cube)[0]
            position = list(position)
            position.insert(0, elasticity)
            writer.writerow(position)       

            time.sleep(1. / 3640.) 

        p.resetBasePositionAndOrientation(cube, [0, 0, 1], p.getQuaternionFromEuler([0, 0, 0]))

    return


def test_orientation(cube):

    fields = ['Theta Y', 'X', 'Y', 'Z', 'Theta X', 'Theta Y', 'Theta Z']

    file = "Tests/test_orientation.csv"
    csvfile = open(file, 'w')
    writer = csv.writer(csvfile)
    writer.writerow(fields)

    intervals = numpy.linspace(0, numpy.pi / 2, 10)

    for i in range(0, 10):

        initial_position = [0.0, 0.0, 1.0]
        theta_y_init = round(intervals[i], 3)
        initial_orienation = [0.0, theta_y_init, 0.0]

        initial_position.insert(0, theta_y_init)
        initial_position.extend(initial_orienation)
        writer.writerow(initial_position)

        p.resetBasePositionAndOrientation(cube, [0, 0, 1], p.getQuaternionFromEuler([0, theta_y_init, 0]))

        for t in range(0, 200):

            p.stepSimulation()

            position, quaternion = p.getBasePositionAndOrientation(cube)
            orientation = p.getEulerFromQuaternion(quaternion)

            position = list(position)
            position.insert(0, theta_y_init)
            position.extend(orientation)
            
            for i in range(len(position)):
                position[i] = round(position[i], 3)

            writer.writerow(position)        

            time.sleep(1. / 3640.)

    return
