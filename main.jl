# Galileo 3
#
# Infers the elasticity of a bouncing soft body from a sequence of postition observations and predicts its future trajectory 
# What representation in the generative model produces inferences that correlate best with human judgments?
#      * using a rigid body
#      * using a sphere
#      * etc.
# Ground truth trajectories are generated in RealFlow and read from .csv files
#
# TODO: Parallel implementation
# TODO: Infer elasticity for trajectories simulated in PyBullet using both MCMC and SMC

using Distributed
using ZipFile

include("args.jl")
include("Model/bouncing_object.jl")
include("Inference/mcmc.jl")
include("Inference/particle_filter.jl")


parallel = false
debug = false

@everywhere begin
function run(fname, args, w1, w2)
    
    # Initialize Bullet simulation context 
    debug_viz = false
    if args.algorithm=="DEBUG"
        debug_viz = true
        client = bullet.connect(bullet.GUI)::Int64
        bullet.resetDebugVisualizerCamera(2, 0, 0, [0, 0, .5])
    else 
        client = bullet.connect(bullet.DIRECT)::Int64
    end
    
    bullet.setAdditionalSearchPath(pybullet_data.getDataPath())

    if args.gt_source=="Bullet"
        sim = BulletSim(step_dur=1/60; client=client)
        fname = string("Tests/BulletStimulus/Data/", args.gt_shape, "/", fname) # Good
    else 
        sim = BulletSim(step_dur=1/30; client=client)
        fname = string("Data/RealFlowData/", args.expt_id, "/", fname)          # Good
    end


    # Initialize scene, including target object state (using observation data)
    init_scene(debug_viz)
    init_position, _, init_velocity, observations, t_s = read_obs_file(fname, args.algorithm)
    init_state = init_target_state(sim, args.target_id, init_position, init_velocity)
    sim_args = (sim, init_state, t_s)

    
    # Generate a sample trace (no inference) 
    if args.algorithm == "DEBUG"

        # Generate and display a single trace 
        println("DEBUGGING")
        (trace, _) = Gen.generate(generate_trajectory, sim_args, observations)
        display(Gen.get_choices(trace))

    # Run inference using MCMC
    elseif args.algorithm == "MCMC"
        
        println("Initializing MCMC")
        avg_last_hundred, results = gaussian_drift_inference(generate_trajectory, sim_args, observations)
        println("Elasticity estimate: ", avg_last_hundred)

        # Write output trajectory to a .csv file
        tokens = split(fname, "_")
        tokens = split(fname, "_")
        fname = string(tokens[2], "_", tokens[3])
        println("Writing to: ", fname)
        f = ZipFile.addfile(w1, fname)
        write_to_csv(results, f)
        
    # Run inference using particle filter    
    else

        # Filter n particles to explain the complete trajectory
        println("Initializing Particle Filter")
        n = args.num_particles
        s = args.save_intermediate
        results, _ = infer(generate_trajectory, sim_args, observations, 
                            n, s, w2, fname)

        # Write output particles to a .csv
        tokens = split(fname, "_")
        fname = string(tokens[2], "_", tokens[3])
        f = ZipFile.addfile(w1, fname)
        write_to_csv(results, f)

        if args.predict

            # For each particle, simulate next 90 time steps
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

    # Extract ground truth trajectory filenames from correct directory
    args.gt_source == "RealFlow" ? 
        dir = string("Data/RealFlowData/", args.expt_id, "/") : 
        dir = string("Tests/BulletStimulus/Data/", args.gt_shape, "/")
    fnames = readdir(dir)
    fnames = filter_unwanted_filenames(fnames)
    sort!(fnames, lt=trial_order)

    # Generate output filepath(s) and writers
    noise_id = generate_noise_id(args)
    output_path = string(args.expt_id, "/", args.model_id, "/", args.target_id, "/", noise_id, "/", args.algorithm)
    w1, w2 = make_directories_and_writers(output_path)

    # Execute particle filter on all input trajectories
    if parallel         # Distribute to workers
        if nworkers() == 1
            addprocs(15)
        end
        pmap(fname -> run(fname, args, w1, w2), fnames)
    else
        if debug        # Run one test file
            run(fnames[1], args, w1, w2)
        else
            for fname in fnames
                run(fname, args, w1, w2)
            end
        end
    end

    close(w1)
    if args.algorithm == "SMC"
        close(w2)  
    end
end

main()
