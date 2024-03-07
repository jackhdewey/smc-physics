# For each of a set of files describing 'ground truth' trajectories, attempt to infer the elasticity of the bouncing object and predict the next 60 frames
# Given its elasticity and geometry, how will the object bounce (what will be its trajectory)?
#
# GOAL: Perform sequential inference of elasticity using 30 frames of observations, followed by 60 frames of prediction
# GOAL: Test what simulation settings (i.e. resource-rational approximations) best match human performance:
#           * using a rigid body
#           * using a sphere
#           * adding noise to orientation
#           * etc.
#
# DONE: Print / write the traces to a .csv file, check whether it recovers elasticity
# DONE: Save itermediate particle filter state
#
# TODO: Fix off-by-one error

include("inference.jl")
include("Utilities/plots.jl")

function main()

    # Read ground truth trajectories
    dir = "RealFlowData/"
    fnames = readdir(dir)
    fnames = filter_unwanted_filenames(fnames)
    println(fnames)

    for i in eachindex(fnames)

        # Initialize simulation context 
        client = bullet.connect(bullet.DIRECT)::Int64
        bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
        #bullet.resetSimulation(bullet.RESET_USE_DEFORMABLE_WORLD)
        bullet.resetDebugVisualizerCamera(4, 0, -4, [0, 0, 2])
        sim = BulletSim(step_dur=1/30; client=client)

        # Initialize simulation using observed data
        fname = fnames[i]
        initial_position, initial_velocity, observations = read_observation_file(fname)
        init_state = generate_scene(sim, "cube", initial_position, initial_velocity)
        args = (30, init_state, sim)

        # Simulate 20 particles to estimate elasticity from observed trajectory 
        results, weights = infer(fname, args, observations, 20)
        gif(animate_traces(results), fps=24)

        # Write out put particles to a .csv file
        tokens = split(fname, "_")
        fname = string("BulletData/Observations/observations_", tokens[2], "_", tokens[3])
        write_to_csv(results, fname)
            
        # For each particle, predict the next 90 time steps
        ppd = predict(results, 90)
        gif(animate_traces(ppd), fps=24)

        # Write predicted trajectories to a .csv file
        fname = string("BulletData/Predictions/predictions_", tokens[2], "_", tokens[3])
        write_to_csv(ppd, fname)

        bullet.disconnect()
    end
end


main()