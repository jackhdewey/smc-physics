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

include("Model/bouncing_object.jl")
include("Inference/particle_filter.jl")
include("Utilities/plots.jl")


function main()

    # Select the target object type(s)
    stimulus_id = "Exp2"
    output_id = "Exp2/Cube"
    target_id = "Cube"

    # Number of particles to be sampled during inference
    num_particles = 20

    # Read ground truth trajectories
    dir = string("RealFlowData/", stimulus_id, "/")
    fnames = readdir(dir)
    fnames = filter_unwanted_filenames(fnames)
    sort!(fnames, lt=trial_order)

    #=
    provide file names as a csv / text file
    read them into a dictionary / set
    if (dictionary.contains(names[i])) 
        else 
    =#

    # Possibly add ability to run a specific trial / file
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

        # Filter 20 particles through the complete trajectory
        num_timesteps = length(observations)
        args = (sim, init_state, num_timesteps)
        results, weights = infer(fname, output_id, generate_trajectory, args, observations, num_particles)
        
        # For each output particle, predict the next 90 time steps
        ppd = predict(results, 90)

        #gif(animate_traces(results), fps=24)
        #gif(animate_traces(ppd), fps=24)

        # Write inferred trajectories to a .csv file
        tokens = split(fname, "_")
        fname = string("BulletData/", output_id, "/Inferences/inferences_", tokens[2], "_", tokens[3])
        write_to_csv(results, fname)

        # Write predicted trajectories to a .csv file
        fname = string("BulletData/", output_id, "/Predictions/predictions_", tokens[2], "_", tokens[3])
        write_to_csv(ppd, fname)

        bullet.disconnect()
    end
end


main()
