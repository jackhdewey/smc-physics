import pybullet as p
import pybullet_data
import simulation

# TODO: Find a way to run multiple simulations, save, display as 'video' later
def main():
    p.connect(p.GUI)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())
    p.resetDebugVisualizerCamera(4, 0, -1, [0, 0, 1])

    spheres, locations = simulation.init_scene()
    simulation.test_orientation(spheres[0])

    p.disconnect()
    return

main()