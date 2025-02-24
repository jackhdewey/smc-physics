# Data analysis
# Evaluate and visualize model performance (at inferring elasticity, relative to ground truth and/or human judgments)
#   Generate 2D plots showing model-inferred elasticities vs human judgments
#   Identify trials with poor model / human correlation by absolute error

using DataFrames
using DataFramesMeta
using Statistics
using Plots

include("../args.jl")
include("utils.jl")


#########
# PLOTS #
#########

# Plot model against ground truth
function plot_vs_gt(type, data, expt_id, target_id, marker_shape, plots_path)

    if type == "sim"
        ylabel = "Model Estimate"
        title = "$expt_id Model ($target_id) v. Ground Truth Individual Stimuli"
    else 
        ylabel = "Human Estimate"
        title = "$expt_id Human v. Ground Truth Individual Stimuli"
    end

    # default(aspect_ratio = :equal)
    p = palette(:jet)
    scatter(data.gtElasticity,
            data.judgment,
            #yerror=sim.std_err_mean ./ 2,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor = data.gtElasticity,          # zcolor = :gtElasticity,
            clims=(0, 1),
            xlims=(0, 1.05),
            ylims=(0, 1),
            xlabel = "Ground Truth Elasticity",
            ylabel = ylabel,
            title = title,
            colorbar=true,
            legend=false,
            palette=p,
            markershape = marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(data.gtElasticity, data.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)

    type == "sim" ?
        savefig(joinpath(plots_path, "sim_vs_gt_individual_stimuli.png")) :
        savefig(joinpath(plots_path, "human_vs_gt_individual_stimuli.png"))

end


# Plot mean model estimates against mean human judgments
function plot_sim_vs_human(sim_data, human_data, expt_id, target_id, marker_shape, plots_path)

    yerror = human_data.std_err_mean ./ 2

    title = "$expt_id Model ($target_id) vs. Human Individual Stimuli"
    filename = "sim_vs_human_individual_stimuli"

    # @infiltrate
    # default(aspect_ratio = :equal)
    p = palette(:jet)
    scatter(human_data.judgment,
            sim_data.judgment,
            yerror = yerror,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor = human_data.gtElasticity,
            xlims=(0, 1.05),
            ylims=(0, 1),
            clims=(0, 1),
            xlabel="Human",
            ylabel="Model",
            title=title,
            legend=false,
            colorbar=true,
            palette=p,
            markershape=marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(sim_data.judgment, human_data.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)

    savefig(joinpath(plots_path, filename))


    title = "$expt_id Model ($target_id) vs. Human Mean"
    filename = "sim_vs_human_average_judgments"

    sim_data = @chain sim_data begin
        @groupby :gtElasticity
        @combine :judgment = mean(:judgment)
    end

    human_data = @chain human_data begin
        @groupby :gtElasticity
        @combine :judgment = mean(:judgment)
    end
    # gui()

    # @infiltrate
    # default(aspect_ratio = :equal)
    p = palette(:jet)
    scatter(human_data.judgment,
            sim_data.judgment,
            yerror = yerror,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor = human_data.gtElasticity,
            xlims=(0, 1.05),
            ylims=(0, 1),
            clims=(0, 1),
            xlabel="Human",
            ylabel="Model",
            title=title,
            legend=false,
            colorbar=true,
            palette=p,
            markershape=marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(sim_data.judgment, human_data.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)

    savefig(joinpath(plots_path, filename))

end


function main()

    args = Args()

    # Read simulation data into data frame and store high-error trials
    noise_id = generate_noise_id(args)
    inference_param_id = generate_inference_param_id(args)
    sim_data = read_simulation_data(args.expt_id, args.model_id, args.target_id, noise_id, args.algorithm, inference_param_id, false)
    sim_data, error_trials = process_individual_stimuli_sim(sim_data)
    println("High Error Trials: ", error_trials)

    # Generate output filepath
    plots_path = generate_plot_path(args.expt_id, args.model_id, args.target_id, noise_id, inference_param_id, "TestJudgments")
    println("Plots Path: ", plots_path)

    # Plotting variables
    if args.target_id == "Cube"
        marker_shape = :square
    else
        marker_shape = :circle
    end

    # Plot inferred elasticity against ground truth, display correlation
    plot_vs_gt("sim", sim_data, args.expt_id, args.target_id, marker_shape, plots_path)

    if occursin("Exp", args.expt_id)
        # Read human data into data frame
        human_data = read_subject_data(args.expt_id)
        human_data = process_individual_stimuli_human(human_data)

        # Plot human inferred elasticity against ground truth, display correlation
        plot_vs_gt("human", human_data, args.expt_id, args.target_id, marker_shape, plots_path)
        plot_sim_vs_human(sim_data, human_data, args.expt_id, args.target_id, marker_shape, plots_path)
    end

end

main()
