using CSV
using DataFramesMeta
using Glob
using Printf
using StatsPlots

base = "/Users/maxs/smc-physics/"

data_path = joinpath(base, "debug_bulletxbullet/")

infPath = joinpath(data_path, "inferences")
partPath = joinpath(data_path, "intermediate")

inferences = readdir(infPath)
particles = readdir(partPath)

nPar = 20
num_stim = length(inferences)
elast_chars = [split(f, '_')[1][4:end] for f in inferences]
elasticities = [parse(Float64, c) * 0.1 for c in elast_chars]
variations = [split(f, ['_', '.'])[2][4:end] for f in inferences]


plots = []
for iStim = 1:1
    currEla = Int64(elasticities[iStim] * 10)
    currVar = variations[iStim]
    id = string("Ela", currEla, "_Var", currVar)

    stimPath = @sprintf("%s/%s_*.csv", partPath, id)
    fList = glob(id * "_*.csv", partPath)
    currModelPred = fill(NaN, nPar, length(fList))
    for iFrm = 1:length(fList)
        iFrmPath = @sprintf("%s/%s_%d.csv", partPath, id, iFrm)
        T = CSV.read(iFrmPath, DataFrame)
        for iPar = 1:nPar
            val = T.elasticity[findall(T.particle .== iPar)[end]]
            currModelPred[iPar, iFrm] = val
        end
    end

    p = violin(
        currModelPred,
        legend=false,
        reuse=false,
        alpha=0.5,
        xticks=0:5:length(currModelPred),
        yticks=1:0.1:1,
        primary=false
    )
    # p = plot()
    nFrames = size(currModelPred)[2]
    for i = 1:nFrames
        data = currModelPred[:, i]
        scatter!(i * ones(20), data, color="black", markersize=1)
        # xlims!(1, nFrames)
        yticks!(0:0.1:1)

    end

    currElaDecimal = currEla
    println(currElaDecimal)
    hline!([currElaDecimal, currElaDecimal], color="black", label="GT")

    display(p)
    # push!(plots, p)

end

# for p in plots
#     display(p)
# end
