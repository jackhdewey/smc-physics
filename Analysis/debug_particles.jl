using CSV
using Plots
using DataFramesMeta
using DataFrames
using Glob

include("../Utilities/fileio.jl")
include("../Utilities/plots.jl")
dir = "Data/BulletData/Modelv5/Sphere/PosVar075/Test/Intermediate/"
dir = "Data/BulletData/Modelv5/Cube/PosVar075/Test/Intermediate/"

# particle_files = filter_unwanted_filenames(readdir(dir))
particle_files = glob("*Ela0_Var3*", dir)
output_id = "Modelv5/Cube/PosVar075/Test"
dir = string("/Users/maxs/smc-physics/Data/BulletData/", output_id, "/Intermediate/")

num_timesteps = 20
num_particles = 20


for ela = 0:9
    plt = plot()

    gt_file = "/Users/maxs/smc-physics/Data/RealFlowData/Test/Cube_Ela" * "$ela" * "_Var1_observed.csv"
    ground_truth = CSV.read(gt_file, DataFrame)

    num_timesteps = size(ground_truth)[1]
    println(num_timesteps)

    all_data = zeros(20, num_timesteps)
    for t = 1:num_timesteps

        file = string("/Users/maxs/smc-physics/Data/BulletData/Modelv5/Cube/PosVar075/Test/Intermediate/particles_Ela", "$ela", "_Var1_$t.csv")
        println(file)
        data = CSV.read(string(file), DataFrame)

        # For each particle
        for p = 1:num_particles

            particle = data[data.particle.==p, :]
            elasticity = particle[!, :elasticity][1]
            all_data[p, t] = elasticity
        end
    end
    for i = 1:20
        # plot!(plt, all_data[i, :], legend=false, xticks=1:2:num_timesteps)
        plot!(plt, all_data[i, :], legend=false, ylims = (0, 1))
    end
    # display(plt)
    savefig(string("Analysis/Plots/", "Elasticity_0.$ela", ".png"))
end
# for p = 1:20
#     pdata = data[data.particle.==p, :]
#     x = 1:maximum(data.frame)
#     plt = plot(x, pdata.elasticity, legend=false, ylim=(0, 0.2), size=(1200, 800))
#     push!(plots, plt)

# end
# plot(plots..., layout=(4, 5))
# display(p)
