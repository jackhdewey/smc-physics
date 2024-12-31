# Parallel implemention 

using Distributed

include("Utilities/fileio.jl")

parallel = true
debug = false

function main()

    @everywhere include("main.jl")

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
    output_path = make_directories(output_path)
    println(output_path)
    zip = false

    # Execute particle filter on all input trajectories
    if parallel    

        #@everywhere include("main.jl")

        # Distribute to workers
        if nworkers() == 1
            addprocs(15)
        end
        pmap(fname -> run_inference(fname, args, output_path, zip), fnames)
        println("DONE")

    else

        #include("main.jl")

        if debug        
            # Run one test file
            run_inference(fnames[1], args, output_path, zip)
        else
            for fname in fnames
                println(fname)
                run_inference(fname, args, output_path, zip)
            end
        end
        
    end

end

main()
