# Main function for running pybullet tests defined in simulation.py
# TODO: Find a way to save simulations in order to display as video later

import pybullet as p
import pybullet_data
import simulation

def main():
    p.connect(p.GUI)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())
    p.resetDebugVisualizerCamera(3, 0, -1, [0, 0, 1])

    cubes, locations = simulation.init_scene()
    simulation.test_orientation(cubes[0])

    p.disconnect()
    return

main()
