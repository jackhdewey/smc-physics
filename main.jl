# For each of a set of files describing 'ground truth' trajectories, attempt to infer the elasticity of the bouncing object and predict the next 60 frames
# Given its elasticity and geometry, how will the object bounce (what will be its trajectory)?
#
# GOAL: Infer elasticity from 30 frames of observations, then predict the next 60 frames
# GOAL: Test what simulation settings (i.e. resource-rational approximations) best match human performance:
#          * using a rigid body
#          * using a sphere
#          * adding noise to orientation
#          * etc.
#
# DONE: Print / write the traces to a .csv file, check whether it recovers elasticity
# DONE: Save itermediate particle filter state
#
# TODO: Infer spheres using spheres
# TODO: Analyze average performance at inferring elasticity
#
# CONSIDER: Infer cubes using spheres
# CONSIDER: Fix off-by-one error

include("inference.jl")
include("Utilities/plots.jl")

function main()

    # Read ground truth trajectories
    dir = "RealFlowData/Sphere/"
    fnames = readdir(dir)
    fnames = filter_unwanted_filenames(fnames)
    println(fnames)

    #for i in eachindex(fnames)
    for i in range(16):

        # Initialize simulation context 
        client = bullet.connect(bullet.DIRECT)::Int64
        bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
        bullet.resetDebugVisualizerCamera(4, 0, -4, [0, 0, 2])
        sim = BulletSim(step_dur=1/30; client=client)
        init_scene()

        # Initialize target state using observed data
        fname = fnames[1]
        initial_position, initial_velocity, observations = read_observation_file(fname, i)
        init_state = init_target_state(sim, "sphere", initial_position, initial_velocity)

        # Given an initial state and observed trajectory, filter 20 particles to estimate elasticity
        args = (sim, init_state, 30)
        results, weights = infer(fname, args, observations, 20)
        gif(animate_traces(results), fps=24)

        # For each output particle, predict the next 90 time steps
        ppd = predict(results, 90)
        #gif(animate_traces(ppd), fps=24)

        # Write inferred trajectories to a .csv file
        tokens = split(fname, "_")
        fname = string("BulletData/Sphere/Observations/observations_", tokens[2], "_", tokens[3])
        write_to_csv(results, fname)

        # Write predicted trajectories to a .csv file
        fname = string("BulletData/Sphere/Predictions/predictions_", tokens[2], "_", tokens[3])
        write_to_csv(ppd, fname)

        bullet.disconnect()
    #end
end


main()
