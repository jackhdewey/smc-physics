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

using Distributed
using Parameters
using ZipFile

include("Model/bouncing_object.jl")
include("Inference/particle_filter.jl")
include("Utilities/plots.jl")


@with_kw struct Args

    model_id::String = "Modelv5"
    target_id::String = "Cube"
    noise_id::String = "PosVar05"
    expt_id::String = "Test"
    output_id::String = string(model_id, "/", target_id, "/", noise_id, "/", expt_id)

    debug::Bool = false

    # Inference parameters
    num_particles::Int = 20
    save_intermediate::Bool = true

    # Prediction parameters
    predict::Bool = false
    prediction_timesteps::Int = 90

end

function main()

    args = Args()

    # Read ground truth trajectories
    dir = string("Data/RealFlowData/", args.expt_id, "/")
    fnames = readdir(dir)
    fnames = filter_unwanted_filenames(fnames)
    sort!(fnames, lt=trial_order)

    output_id = string(args.output_id, "1")
    dir_base = string("Data/BulletData/", output_id)
    if !isdir(dir_base)
        mkdir(dir_base)
    end

    particle_dir = string(dir_base, "/Intermediate/")
    if !isdir(particle_dir)
        mkdir(particle_dir)
    end
    w2 = ZipFile.Writer(string(particle_dir, "/particles.zip"))

    output_dir = string(dir_base, "/Inferences/")
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    w1 = ZipFile.Writer(string(output_dir, "/inferences.zip"))
    
    for fname in fnames
        run(fname, args, w1, w2)
    end

    # Map each observed trajectory to a process executing a particle filter
    #pmap(fname -> run(fname, args), fnames)

    close(w1)
    close(w2)
        
end

@everywhere begin
function run(fname, args, w1, w2)
    
    # Initialize simulation context 
    client = bullet.connect(bullet.DIRECT)::Int64
    bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
    bullet.resetDebugVisualizerCamera(4, 0, -4, [0, 0, 2])
    sim = BulletSim(step_dur=1/30; client=client)

    # Initialize target object state using observed data
    init_scene()
    fname = string(args.expt_id, "/", fname)
    init_position, _, init_velocity, observations, t_s = read_obs_file(fname, args.debug)
    init_state = init_target_state(sim, args.target_id, init_position, init_velocity)
    model_args = (sim, init_state, t_s)

    if args.debug
        (trace, _) = Gen.generate(generate_trajectory, model_args, observations)
        display(Gen.get_choices(trace))
    else

        # Filter particles to explain the complete trajectory
        results, _ = infer(generate_trajectory, model_args, observations, w2, args.num_particles, args.save_intermediate, fname)

        # Write output particles to a .csv file
        tokens = split(fname, "_")
        fname = string(tokens[2], "_", tokens[3])
        f = ZipFile.addfile(w1, fname)
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

    #=
    catch e
        println(e) 
        println("Disconnecting Bullet")
        bullet.disconnect()
    end
    =#

end
end

main()
