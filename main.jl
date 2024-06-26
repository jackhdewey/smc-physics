# Galileo 3
#
# Infers the elasticity of a bouncing object from n observation frames and predicts its future trajectory 
# Given how it has bounced so far, what will be the target object's future trajectory?
#
# What kind of generative model correlates best with human judgments?
#       * using a rigid body
#       * using a sphere
#       * etc.
#
# Ground truth trajectories are generated in RealFlow and read from .csv files

include("Model/bouncing_object.jl")
include("Inference/particle_filter.jl")
include("Utilities/plots.jl")


function main()

    # Select the target object type(s)
    model_id = "Modelv5/PosVar075"
    expt_id = "Exp1"
    target_id = "Sphere"
    output_id = string(model_id, "/", expt_id, "/", target_id)

    # Inference parameters
    num_particles = 20
    save_intermediate = true

    # Prediction parameters
    prediction_timesteps = 90

    # Read ground truth trajectories
    dir = string("RealFlowData/", expt_id, "/")
    fnames = readdir(dir)
    fnames = filter_unwanted_filenames(fnames)
    sort!(fnames, lt=trial_order)

    # Iterate through each observed trajectory, executing a corresponding particle filter to infer elasticity 
    for i in eachindex(fnames)
        
        # Initialize simulation context 
        client = bullet.connect(bullet.DIRECT)::Int64
        bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
        bullet.resetDebugVisualizerCamera(4, 0, -4, [0, 0, 2])
        sim = BulletSim(step_dur=1/30; client=client)

        # Initialize target object state using observed data
        init_scene()
        fname = string(expt_id, "/", fnames[i])
        initial_position, initial_orientation, initial_velocity, observations = read_observation_file(fname)
        init_state = init_target_state(sim, target_id, initial_position, initial_velocity)

        # Filter particles through the complete trajectory
        num_timesteps = length(observations)
        args = (sim, init_state, num_timesteps)
        
        # Filter particles to explain the complete trajectory
        results, _ = infer(fname, output_id, generate_trajectory, args, observations, num_particles, save_intermediate)
    
        # Write inferred trajectories to a .csv file
        tokens = split(fname, "_")
        fname = string("BulletData/", output_id, "/Inferences/inferences_", tokens[2], "_", tokens[3])
        write_to_csv(results, fname)

        #=
        # For each output particle, predict the next 90 time steps
        ppd = predict(results, prediction_timesteps)
        #gif(animate_traces(ppd), fps=24)

        # Write predicted trajectories to a .csv file
        fname = string("BulletData/", output_id, "/Predictions/predictions_", tokens[2], "_", tokens[3])
        write_to_csv(ppd, fname)
        =#

        bullet.disconnect()

    end
end


main()
