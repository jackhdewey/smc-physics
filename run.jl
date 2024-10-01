# Parallel implementation

using Distributed
if nworkers() == 1
    addprocs(15)
end

include("Utilities/fileio.jl")

@everywhere include("main.jl")
# include("main.jl")

all_args = []
all_fnames = []

for i = 2:4
    expt_id = "Exp$i"
    gt_source = "RealFlow"
    args = Args(expt_id=expt_id, gt_source=gt_source, target_id="Sphere", save_intermediate=true)

    if gt_source == "RealFlow"
        dir = string("Data/RealFlowData/", args.expt_id, "/")
    else
        dir = string("Tests/BulletStimulus/Data/")
    end

    expt_fnames = readdir(dir)
    expt_fnames = filter_unwanted_filenames(expt_fnames)

    append!(all_fnames, expt_fnames)
    append!(all_args, fill(args, length(expt_fnames)))

end

# sort!(all_fnames, lt=trial_order)
pmap(run, zip(all_fnames, all_args))
