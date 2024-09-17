# Galileo 3
#
# Infers the elasticity of a bouncing soft body from a sequence of postition observations and predicts its future trajectory 
# What representation in the generative model produces inferences that correlate best with human judgments?
#      * using a rigid body
#      * using a sphere
#      * etc.
# Ground truth trajectories are generated in RealFlow and read from .csv files
#
# TODO: Try inference using MCMC
# TODO: Infer elasticity for trajectories simulated in PyBullet
# TODO: Display particle trajectories perspicuously 

using Distributed
using Parameters
using ZipFile

include("Model/bouncing_object.jl")
include("Inference/mcmc.jl")
include("Inference/particle_filter.jl")
include("Utilities/plots.jl")
include("Utilities/args.jl")

@everywhere begin
function run(fname, args, w1, w2)
    
    # Initialize simulation context 
    client = bullet.connect(bullet.DIRECT)::Int64
    bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
    bullet.resetDebugVisualizerCamera(4, 0, -4, [0, 0, 2])
    sim = BulletSim(step_dur=1/30; client=client)

    # Initialize scene 
    init_scene()

    # Initialize target object state using observation data
    args.gt_source == "RealFlow" ? fname = string("Data/RealFlowData/", args.expt_id, "/", fname) : fname = string("Tests/Data/", fname)
    init_position, _, init_velocity, observations, t_s = read_obs_file(fname, args.algorithm == "PARTICLE_FILTER")
    init_state = init_target_state(sim, args.target_id, init_position, init_velocity)
    model_args = (sim, init_state, t_s)

    if args.algorithm == "DEBUG"

        # Generate and display a single trace 
        println("DEBUGGING")
        (trace, _) = Gen.generate(generate_trajectory, model_args, observations)
        display(Gen.get_choices(trace))

    elseif args.algorithm == "MCMC"
        
        println("Initializing MCMC")
        trace = gaussian_drift_inference(model_args, observations)
        println(trace[:latents => 1 => :restitution])

        # Write output trajectory to a .csv file
        tokens = split(fname, "_")
        tokens = split(fname, "_")
        fname = string(tokens[2], "_", tokens[3])
        println(fname)
        #f = ZipFile.addfile(w1, fname)
        #write_to_csv(results, f)
        
    else

        println("Initializing Particle Filter")

        # Filter n particles to explain the complete trajectory
        results, _ = infer(generate_trajectory, model_args, observations, w2, args.num_particles, args.save_intermediate, fname)

        # Write output particles to a .csv file
        tokens = split(fname, "_")
        fname = string(tokens[2], "_", tokens[3])
        println(fname)
        f = ZipFile.addfile(w1, fname)
        write_to_csv(results, f)

        if args.predict

            # For each particle, simulate the next 90 time steps
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

########
# MAIN #
########

function main()

    args = Args()

    # Read ground truth trajectories
    args.gt_source == "RealFlow" ? dir = string("Data/RealFlowData/", args.expt_id, "/") : dir = string("Tests/BulletStimulus/Data/")
    fnames = readdir(dir)
    fnames = filter_unwanted_filenames(fnames)
    sort!(fnames, lt=trial_order)

    w1, w2 = make_directories_and_writers(args)

    for fname in fnames
        run(fname, args, w1, w2)
    end

    # Map each observed trajectory to a process executing a particle filter
    #pmap(fname -> run(fname, args, w1, w2), fnames)

    close(w1)
    close(w2)  
end

main()
