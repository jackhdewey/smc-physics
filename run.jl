# Parallel implemention

using Distributed
if nworkers() == 1
    addprocs(3)
end

include("Utilities/fileio.jl")

@everywhere include("main.jl")
# include("main.jl")

all_args = []
for i = 2:4
    expt_id = "Exp$i"
    args = Args(expt_id=expt_id, gt_source="RealFlow", target_id="Sphere")

    # append!(all_fnames, [(joinpath(expt_id, f), args) for f in fnames])
    push!(all_args, args)

end


# println(all_args)
pmap(main, all_args)
# map(main, all_args)
