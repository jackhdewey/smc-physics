# Main function for running PyBullet tests defined in simulation.py

# TODO: Find a way to save simulations in order to display as video later

import pybullet as p
import pybullet_data
import tests

def main():
    p.connect(p.GUI)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())
    p.resetDebugVisualizerCamera(2, 0, -1, [0, 0, .5])

    rigidBodies = tests.init_scene()
    tests.test_orientation(rigidBodies[0])

    p.disconnect()
    return

main()
