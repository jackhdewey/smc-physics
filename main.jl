# Galileo 3
#
# Infers the elasticity of a bouncing soft body from a sequence of position observations and predicts its future trajectory
# What representation in the generative model produces inferences that correlate best with human judgments?
#      * using a rigid body
#      * using a sphere
#      * etc.
# Ground truth trajectories are generated in RealFlow and read from .csv files
#
# TODO: Try inference using MCMC
# TODO: Infer elasticity for trajectories simulated in PyBullet
# TODO: Display particle trajectory data

using ZipFile
using Random

# Seed the random number generator
Random.seed!(1234)
include("Model/bouncing_object.jl")
include("Inference/mcmc.jl")
include("Inference/particle_filter.jl")
include("Utilities/plots.jl")
include("Utilities/args.jl")

function run(fname, args)

    # Initialize simulation context
    client = bullet.connect(bullet.DIRECT)::Int64
    bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
    bullet.resetDebugVisualizerCamera(4, 0, -4, [0, 0, 2])
    sim = BulletSim(step_dur=1 / 30; client=client)

    # Initialize scene 
    init_scene()

    # Initialize target object state using observed data
    fname = string(args.expt_id, "/", fname)
    init_position = nothing
    init_velocity = nothing
    observations = nothing
    t_s = nothing

    try
        init_position, _, init_velocity, observations, t_s = read_obs_file(fname, args, true)
    catch e
        println("---------------------------------------------------------\n\n")
        println(e)
        println("---------------------------------------------------------\n\n")
        return
    end

    init_state = init_target_state(sim, args.target_id, init_position, init_velocity)
    model_args = (sim, init_state, t_s)

    if args.debug
        # Generate and display a single trace
        (trace, _) = Gen.generate(generate_trajectory, model_args, observations)
        display(Gen.get_choices(trace))

    else
        # Filter n particles to explain the complete trajectory
        println("Initializing Particle Filter")
        results, _ = infer(generate_trajectory, model_args, observations, args, args.num_particles, args.save_intermediate, fname)

        # Write output particles to a .csv file
        # println(fname)
        tokens = split(fname, "_")
        fname = string(tokens[2], "_", tokens[3])
        # f = ZipFile.addfile(w1, fname)

        tokens = split(fname, "_")
        f = string("Data/BulletData/", args.output_id, "/Inferences/", tokens[1], "_", tokens[2])
        write_to_csv(results, f)

        if args.predict
            # For each output particle, predict the next 90 time steps
            ppd = predict(results, args.prediction_timesteps)
            #gif(animate_traces(ppd), fps=24)

            # Write predicted trajectories to a .csv file
            fname = string(dir_base, "/Predictions/predictions_", tokens[2], "_", tokens[3])
            write_to_csv(ppd, fname)
        end
    end

    bullet.disconnect()
end
