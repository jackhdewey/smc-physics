# Starting point for inference using a probabilistic generative model, with options to:
#    1. debug the model (run and display a single trace)
#    2. run inference in parallel or on a single process

using Distributed
if nworkers() == 1
    addprocs(15)
end
@everywhere include("Inference/run_inference.jl")

using ZipFile
include("Utilities/fileio.jl")

debug = false
parallel = true
zip = false

function main()

    # Instance of modeling and inference arguments
    args = Args()

    # Extract ground truth trajectory files from relevant directory
    if contains(args.expt_id, "Bullet")  
        bullet_shape = split(args.expt_id, "_")[2] 
        dir = string("Data/BulletStimulus/New/", bullet_shape, "/") 
        fnames = readdir(dir)
    else
        dir = string("Data/RealFlowData/", args.expt_id, "/") 
        fnames = filter_unwanted_filenames(readdir(dir))
    end
    sort!(fnames, lt=trial_order)

    # Generate output filepath(s)
    noise_id = generate_noise_id(args)
    inference_id = generate_inference_param_id(args)
    output_path = make_directories(joinpath(args.expt_id, args.model_id, args.target_id, noise_id, args.algorithm, inference_id, "Data"))
    println("Output Filepath... ", output_path)

    # Generate writers (not currently used)
    #w1, w2 = make_writers(output_path, args.algorithm)
    w1, w2 = nothing, nothing

    # Iterate over all parameter settings
    for o_noise in args.observation_noise
        for t_noise in args.transition_noise
            for n_p in args.num_particles
                for r_m in args.rejuvenation_moves
                    params = [o_noise, t_noise, n_p, r_m]

                    # Execute particle filter on all observed trajectories
                    if parallel    
                        # Distribute to workers
                        pmap(fname -> run_inference(fname, args, output_path, w1, w2), fnames)
                        println("DONE")
                    else

                        if debug                 
                            # Run one test file
                            run_inference(fnames[1], args, output_path, w1, w2)
                        else                     
                            # Run all files on one process
                            for fname in fnames
                                println(fname)
                                run_inference(fname, args, output_path, w1, w2)
                            end
                        end
                        
                    end

                end
            end
        end
    end
    
  
    #=
    # Close writers
    close(w1)
    if args.algorithm == "SMC"
        close(w2)  
    end
    =#

end

main()
