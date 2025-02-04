# Starting point for performing inference with a generative model in the loop, with options to:
#    1. run in parallel or on a single process
#    2. debug (run and display a single trace of the generative model)

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

    args = Args()

    # Extract ground truth trajectory files from appropriate directory
    if contains(args.expt_id, "Bullet")  
        bullet_shape = split(args.expt_id, "_")[2] 
        dir = string("Tests/BulletStimulus/Data/", bullet_shape, "/") 
    else
        dir = string("Data/RealFlowData/", args.expt_id, "/") 
    end
    fnames = filter_unwanted_filenames(readdir(dir))
    sort!(fnames, lt=trial_order)

    # Generate output filepath(s)
    noise_id = generate_noise_id(args)
    inference_id = generate_inference_param_id(args)
    output_path = make_directories(joinpath(args.expt_id, args.model_id, args.target_id, noise_id, args.algorithm, inference_id, "Data"))
    println("Output Filepath... ", output_path)

    # Generate writers (not currently used)
    #w1, w2 = make_writers(output_path, args.algorithm)
    w1, w2 = nothing, nothing

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
  
    #=
    # Close writers
    close(w1)
    if args.algorithm == "SMC"
        close(w2)  
    end
    =#

end

main()
