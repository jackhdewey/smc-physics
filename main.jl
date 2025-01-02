# Parallel implemention 

using Distributed
using ZipFile

include("Utilities/fileio.jl")

@everywhere include("Inference/run_inference.jl")

debug = false
parallel = true
zip = false

function main()

    @everywhere include("Inference/run_inference.jl")

    args = Args()

    # Extract ground truth trajectory filenames from correct directory
    args.gt_source == "RealFlow" ? 
        dir = string("Data/RealFlowData/", args.expt_id, "/") : 
        dir = string("Tests/BulletStimulus/Data/", args.gt_shape, "/")
    fnames = readdir(dir)
    fnames = filter_unwanted_filenames(fnames)
    sort!(fnames, lt=trial_order)

    # Generate output filepath(s)
    noise_id = generate_noise_id(args)
    output_path = string(args.expt_id, "/", args.model_id, "/", args.target_id, "/", noise_id, "/", args.algorithm)
    output_path = make_directories(output_path)
    w1, w2 = make_writers(output_path, args.algorithm)

    # Execute particle filter on all input trajectories
    if parallel    

        # Distribute to workers
        if nworkers() == 1
            addprocs(15)
        end

        pmap(fname -> run_inference(fname, args, output_path, w1, w2), fnames)
        println("DONE")

    else

        if debug                 # Run one test file
            run_inference(fnames[1], args, output_path, w1, w2)
        else                     # Run all files on one process
            for fname in fnames
                println(fname)
                run_inference(fname, args, output_path, w1, w2)
            end
        end
        
    end
  
    close(w1)
    if args.algorithm == "SMC"
        close(w2)  
    end

    # COMPRESS OUTPUT FILES

end

main()
