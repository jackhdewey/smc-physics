using CSV
using FileIO, ImageMagick
using Statistics
using Plots
using StatsPlots
using DataFramesMeta
using DataFrames

include("../Utilities/fileio.jl")
include("../Utilities/plots.jl")
include("utils.jl")

project_path = dirname(Base.active_project())

dir = joinpath(project_path, "Data/BulletData/Modelv5/Sphere/PosVar075/Test/Intermediate/")
dir = joinpath(project_path, "Data/BulletData/Modelv5/Cube/PosVar075/Test/Intermediate/")

# particle_files = filter_unwanted_filenames(readdir(dir))
output_id = "Modelv5/Cube/PosVar05/Exp1"
dir = joinpath(project_path, "Data/BulletData", output_id, "Intermediate")

num_particles = 20

heuristics = CSV.read(joinpath(project_path, "Analysis", "heuristicPrediction_Exp1.csv"), DataFrame)
humans = read_subject_data(1)

# global all_data
variation = 1
# for ela = 0:9
ela = 0
for variation = 1:10
    scale = 1
    dpi = 1000
    plt = Plots.plot(size=(800 * scale, 600 * scale), dpi=dpi)

    gt_file = joinpath(project_path, "Data/RealFlowData/Exp1/Cube_Ela$ela" * "_Var$variation" * "_observed.csv")
    println(gt_file)
    ground_truth = CSV.read(gt_file, DataFrame)

    num_timesteps = size(ground_truth)[1]

    global all_data = zeros(20, num_timesteps)
    xs = 1:num_timesteps
    for t = 1:num_timesteps

<<<<<< Updated upstream
        file = string("Data/BulletData/Modelv5/Cube/PosVar075/Test/Intermediate/particles_Ela", "$ela", "_Var1_$t.csv")
        println(file)
        =======
        file = string(project_path, "/Data/BulletData/Modelv5/Cube/PosVar05/Exp1/Intermediate/Ela", "$ela", "_Var$variation", "_$t.csv")
>>>>>>> Stashed changes
        data = CSV.read(string(file), DataFrame)

        # For each particle
        for p = 1:num_particles

            particle = data[data.particle.==p, :]
            elasticity = particle[!, :elasticity][1]
            all_data[p, t] = elasticity
        end
    end

    scenario = string("Cube_Ela", ela, "_Var$variation")
    heuristic = heuristics[findall(heuristics[!, "ID"] .== scenario), "heuristic"][1] + 0.5
    human = @chain humans begin
        @subset :elasticity .== ela / 10
        @combine begin
            :judgment = mean(:rating)
        end
    end
    human_judgment = first(human.judgment)

    for i = 1:20
        # plot!(plt, all_data[i, :], legend=false, xticks=1:2:num_timesteps)
        plot!(
            plt,
            xs,
            all_data[i, :],
            ylims=(0, 1),
            xlabel="Time Step",
            ylabel="Estimated Elasticity",
            label="",
            legend=true
        )
    end
    hline!(plt,
        heuristic .* ones(num_timesteps),
        linestyle=:dash,
        color=:red,
        label="heuristic"
    )
    hline!(plt,
        human_judgment .* ones(num_timesteps),
        linestyle=:dash,
        color=:orange,
        legend=true,
        label="human"
    )
    hline!(plt,
        0.1 * ela .* ones(num_timesteps),
        linestyle=:dash,
        color=:purple,
        legend=true,
        label="ground truth"
    )
    pdfname = string("Analysis/Plots/", "Elasticity_0.$ela", "_Var$variation", ".pdf")
    pdfname = joinpath(project_path, pdfname)
    println("pdf name", pdfname)
    Plots.savefig(pdfname)
    println("saved pdf")

    dpi = 2000

    jpgname = string("Analysis/Plots/", "Elasticity_0.$ela", "_Var$variation", ".jpg")
    println(jpgname)
    @async run(`convert -density $dpi -quality 100 $pdfname $jpgname`)

    plt = Plots.plot(size=(800 * scale, 600 * scale), dpi=dpi)

    # scenario = string("Cube_Ela", ela, "_Var$variation")
    # heuristic = heuristics[findall(heuristics[!, "ID"] .== scenario), "heuristic"][1] + 0.5
    # human = @chain humans begin
    #     @subset :elasticity .== ela / 10
    #     @combine begin
    #         :judgment = mean(:elasticity)
    #     end
    # end
    # human = first(human.judgment)
    for i = 1:num_timesteps
        violin!(
            plt,
            xs,
            all_data[:, i],
            legend=false,
            ylims=(0, 1),
            xlabel="Time Step",
            ylabel="Estimated Elasticity",
            label="",
            color=:teal
        )
    end
    hline!(plt,
        heuristic .* ones(num_timesteps),
        linestyle=:dash,
        color=:red,
        label="heuristic"
    )
    hline!(plt,
        human_judgment .* ones(num_timesteps),
        linestyle=:dash,
        color=:orange,
        legend=true,
        label="human"
    )
    hline!(plt,
        0.1 * ela .* ones(num_timesteps),
        linestyle=:dash,
        color=:purple,
        legend=true,
        label="ground truth"
    )
    # add heuristic line


    pdfname = string("Analysis/Plots/", "Elasticity_violin_0.$ela", "_Var$variation", ".pdf")
    Plots.savefig(pdfname)

    jpgname = string("Analysis/Plots/", "Elasticity_violin_0.$ela", "_Var$variation", ".jpg")
    @async run(`convert -density $dpi -quality 100 $pdfname $jpgname`)
end

plots_path = joinpath(project_path, "Analysis", "Plots")
run(`sh -c "rm $plots_path/*.pdf"`)

# for p = 1:20
#     pdata = data[data.particle.==p, :]
#     x = 1:maximum(data.frame)
#     plt = plot(x, pdata.elasticity, legend=false, ylim=(0, 0.2), size=(1200, 800))
#     push!(plots, plt)

# end
# plot(plots..., layout=(4, 5))
# display(p)
