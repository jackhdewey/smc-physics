# Experiment - Galileo 3
#
# Infers the elasticity and predicts the future trajectory of a bouncing object from 1 second (30 frames) of observation
# 'Ground truth' trajectories generated in RealFlow and read from .csv files
#
# QUESTION: Given how it has bounced so far, what will be the target object's future trajectory?
#
# What kind of generative model correlates best with human judgments?
#       * using a rigid body
#       * using a sphere
#       * adding noise to orientation
#       * etc.
#
# TODO: Fix off-by-one error

include("Inference/particle_filter.jl")
include("Utilities/plots.jl")


function main()

    # Select the target object types
    stimulus_id = "Cube"
    target_id = "Sphere"
    output_id = "SpherexCube"

    # Select particle filter parameters
    num_timesteps = 30
    num_particles = 20

    # Read ground truth trajectories
    dir = string("RealFlowData/", stimulus_id, "/")
    fnames = readdir(dir)
    fnames = filter_unwanted_filenames(fnames)
    sort!(fnames, lt=trial_order)

    for i in eachindex(fnames)

        # Initialize simulation context 
        client = bullet.connect(bullet.DIRECT)::Int64
        bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
        bullet.resetDebugVisualizerCamera(4, 0, -4, [0, 0, 2])
        sim = BulletSim(step_dur=1/30; client=client)
        init_scene()

        # Initialize target state using observed data
        fname = string(stimulus_id, "/", fnames[i])
        initial_position, initial_velocity, observations = read_observation_file(fname)
        init_state = init_target_state(sim, target_id, initial_position, initial_velocity)

        # Given an initial state and observed trajectory, filter 20 particles to estimate elasticity
        args = (sim, init_state, num_timesteps)
        results, weights = infer(fname, output_id, args, observations, num_particles)

        # For each output particle, predict the next 90 time steps
        ppd = predict(results, 90)

        #gif(animate_traces(results), fps=24)
        #gif(animate_traces(ppd), fps=24)

        # Write inferred trajectories to a .csv file
        tokens = split(fname, "_")
        fname = string("BulletData/", output_id, "/Observations/observations_", tokens[2], "_", tokens[3])
        write_to_csv(results, fname)

        # Write predicted trajectories to a .csv file
        fname = string("BulletData/", output_id, "/Predictions/predictions_", tokens[2], "_", tokens[3])
        write_to_csv(ppd, fname)

        bullet.disconnect()
    end
end


main()
