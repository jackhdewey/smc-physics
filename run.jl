# Parallel implemention 

using Distributed
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
