# Starting point for performing inference with a generative model in the loop, with options to:
#    1. run in parallel or on a single process
#    2. debug (run and display a single trace of the generative model)

using Distributed

if nworkers() == 1
    # addprocs(15)
    addprocs(4)
end

@everywhere begin
    global args = nothing
    include("Inference/run_inference.jl")
end

using ZipFile
include("Utilities/fileio.jl")


debug = false
parallel = true
zip = false

function main()

    args = Args()

    # Extract ground truth trajectory files from appropriate directory
    contains(args.gt_source, "RealFlow") ? 
        dir = string("Data/RealFlowData/", args.expt_id, "/") : 
        dir = string("Tests/BulletStimulus/Data/", args.gt_shape, "/")
    fnames = filter_unwanted_filenames(readdir(dir))
    sort!(fnames, lt=trial_order)

    # Generate output filepath(s)
    noise_id = generate_noise_id(args)
    inference_id = generate_inference_param_id(args)
    output_path = make_directories(joinpath(args.expt_id, args.model_id, args.target_id, noise_id, args.algorithm, inference_id, "Data"))
    println("Output Filepath... ", output_path)

    # Generate writers (not currently used so null)
    w1, w2 = nothing, nothing
    #w1, w2 = make_writers(output_path, args.algorithm)

    # Execute particle filter on all observed trajectories
    if parallel    
        # Distribute to workers

        arglist = [Args() for i in 1:10]
        for (i, a) in enumerate(arglist)
            a.num_particles = i
        end

        function run_inference_with_process_specific_args(fname, job_args, output_path, w1, w2)
            global args = job_args
            run_inference(fname, output_path, w1, w2)
        end

        function run_inference_all_trials(args)

            noise_id = generate_noise_id(args)
            inference_id = generate_inference_param_id(args)
            output_path = make_directories(joinpath(
                args.expt_id,
                args.model_id,
                args.target_id,
                noise_id,
                args.algorithm,
                inference_id,
                "Data"
            ))
            
            pmap(fname -> run_inference_with_process_specific_args(fname, args, output_path, w1, w2), fnames)
        end

        for a in arglist
            run_inference_all_trials(a)
        end


        println("DONE")
    else

        if debug                 
            # Run one test file
            run_inference(fnames[1], output_path, w1, w2)
        else                     
            # Run all files on one process
            for fname in fnames
                println(fname)
                run_inference(fname, output_path, w1, w2)
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
