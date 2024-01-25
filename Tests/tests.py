import pybullet as p
import pybullet_data
import simulation

# TODO: Find a way to run multiple simulations, save, display as 'video' later
def main():
    p.connect(p.GUI)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())
    p.resetSimulation(p.RESET_USE_DEFORMABLE_WORLD)
    p.resetDebugVisualizerCamera(4, 0, -1, [0, 0, 1])

    sphereBodies, sphere_locations = simulation.init_scene()
    simulation.run_sample(sphereBodies, sphere_locations)

    p.disconnect()
    return

main()