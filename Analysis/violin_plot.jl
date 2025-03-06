# Data analysis
# Generates a 2D plot showing the evolution of elastities in the particle population over time

using CSV
using DataFramesMeta
using Glob
using Printf
using StatsPlots

include("../args.jl")
include("utils.jl")


# Locate relevant data files
base = pwd()
args = Args()
param_id = generate_inference_param_id(args)
data_path = joinpath(base, args.expt_id, args.model_id, param_id, "Data")

inf_path = joinpath(data_path, "inferences")
inferences = readdir(inf_path)

part_path = joinpath(data_path, "intermediate")
particles = readdir(part_path)

# 
elast_chars = [split(f, '_')[1][4:end] for f in inferences]
elasticities = [parse(Float64, c) * 0.1 for c in elast_chars]
variations = [split(f, ['_', '.'])[2][4:end] for f in inferences]

# Generate plots
plots = []
num_particles = args.num_particles
for i = 1:1
    
    ela = Int64(elasticities[i] * 10)
    var = variations[i]
    id = string("Ela", ela, "_Var", var)
    
    files = glob(id * "_*.csv", part_path)
    model_pred = fill(NaN, num_particles, length(files))
    for t_s = 1:length(files)
        var_timestep = @sprintf("%s/%s_%d.csv", part_path, id, t_s)
        particle_data = CSV.read(var_timestep, DataFrame)
        for p = 1:num_particles
            pred_ela = particle_data.elasticity[findall(particle_data.particle .== p)[end]]
            model_pred[p, t_s] = pred_ela
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
    
    nFrames = size(model_pred)[2]
    for i = 1:nFrames
        data = model_pred[:, i]
        scatter!(i * ones(20), data, color="black", markersize=1)
        # xlims!(1, nFrames)
        yticks!(0:0.1:1)
    end

    currElaDecimal = ela
    println(currElaDecimal)
    hline!([currElaDecimal, currElaDecimal], color="black", label="GT")

    display(p)
    # push!(plots, p)

end

# for p in plots
#     display(p)
# end
