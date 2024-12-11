# Parallel implemention 

using Distributed

@everywhere include("main.jl")

include("Utilities/fileio.jl")

parallel = true
debug = false

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

    # Execute particle filter on all input trajectories
    if parallel         

        # Distribute to workers
        if nworkers() == 1
            addprocs(11)
        end
        pmap(fname -> run_inference(fname, args, output_path), fnames)

    else

        if debug        
            # Run one test file
            run_inference(fnames[1], args, output_path)
        else
            for fname in fnames
                run_inference(fname, args, output_path)
            end
        end
        
    end

end

main()

#=
addprocs(4)

include("Utilities/fileio.jl")

@everywhere include("main.jl")

expt_id = "Exp1"
dir = string("RealFlowData/", expt_id, "/")

fnames = readdir(dir)
fnames = filter_unwanted_filenames(fnames)

sort!(fnames, lt=trial_order)

fnames = [joinpath(expt_id, f) for f in fnames]

pmap(main, fnames)
=#
