# Data analysis
# Assess model performance at inferring elasticity relative to ground truth and/or human judgments
# Identify trials with poor model / human correlation by absolute error
# Generate 2D plots showing model-inferred elasticities vs human judgments
#
# TODO: Integrate particle trajectory data 
# TODO: Display MCMC trajectory data

using DataFrames
using DataFramesMeta
using Statistics
using Plots

include("../args.jl")
include("utils.jl")


#################
# Data Analysis #
#################

# Compute error on each stimulus, extract high error trials
function process_individual_stimuli_sim(sim_data)

    # Define model's elasticity judgment as average of output particles
    sim_data_pred = @chain sim_data begin
        @groupby :stimulusID
        @combine begin

            :judgment = mean(:elasticity)       
            :gtElasticity = first(:gtElasticity)  # :elasticity = :gtElasticity

            # compute error w.r.t. ground truth
            :error = mean(:elasticity) - first(:gtElasticity)
        end
        # @subset :elasticity .> 0.6
        @orderby :stimulusID
    end

    # Filter by error threshold, save high error trials
    high_error = sim_data_pred[sim_data_pred.error .> .2, :]
    high_error_trials = high_error[:, 1]

    return sim_data_pred, high_error_trials
end

# Group subject's elasticity judgments by stimulus and compute mean rating for each stimulus
function process_individual_stimuli_human(sub_data)

    nsubs = length(unique(sub_data.filename))   # number of subjects
    sub_data_pred = @chain sub_data begin
        @groupby :stimulusID
        # @groupby :elasticity# :filename
        # @DataFramesMeta.transform :prediction = mean(:rating)     # Gen also has a transform macro
        @combine begin

            :judgment = mean(:rating)
            :std_err_mean = std(:rating) / sqrt(nsubs) 

            :gtElasticity = mean(:elasticity)

        end
        # @subset :gtElasticity .> 0.6
        @orderby :stimulusID
    end

    return sub_data_pred
end


#########
# PLOTS #
#########

# Plot model against ground truth
function plot_sim_vs_gt(sim, marker_shape, plots_path, expt_id, target_id)

    # default(aspect_ratio = :equal)
    p = palette(:jet)
    scatter(sim.gtElasticity,
            sim.judgment,
            #yerror=sim.std_err_mean ./ 2,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor=sim.gtElasticity,          # zcolor = :gtElasticity,
            clims=(0, 1),
            xlims=(0, 1.05),
            ylims=(0, 1),
            xlabel="Ground Truth Elasticity",
            ylabel="Model Estimate",
            title="Exp $expt_id stimulus-specific elasticity ratings - $target_id",
            colorbar=true,
            legend=false,
            palette=p,
            markershape=marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(sim.gtElasticity, sim.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)

    savefig(joinpath(plots_path, "individual_stimuli_model_v_gt.png"))

end

# Plot human judgments against the ground truth
function plot_human_vs_gt(human_data)

    # default(aspect_ratio = :equal)
    p = palette(:jet)
    scatter(human_data.gtElasticity,
        human_data.judgment,
        yerror=human_data.std_err_mean ./ 2,
        aspect_ratio=:equal,
        markersize=5,
        markeralpha=0.5,
        zcolor=human_data.gtElasticity,# zcolor = :gtElasticity,
        clims=(0, 1),
        xlims=(0, 1.05),
        ylims=(0, 1),
        xlabel="Ground Truth Elasticity",
        ylabel="Human Estimate",
        title="Exp $expt_id stimulus-specific human ratings",
        colorbar=true,
        legend=false,
        palette=p,
        markershape=marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(human_data.gtElasticity, human_data.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)
    plots_path = generate_plot_path(expt_id, model_id, target_id, noise_id)
    savefig(joinpath(plots_path, string("individual_stimuli_judgments_", "against_gt_human", "Exp", expt_id, ".png")))
    # gui()

end


# Plot mean model estimates against mean human judgments
function plot_sim_vs_human(sim_data, human_data, plot_mean)

    # default(aspect_ratio = :equal)

    title = "Exp $expt_id stimulus-specific elasticity ratings - $target_id"
    filename = "individual_stimuli_judgments_high_elasticity"


    if plot_mean

        title = "Exp $expt_id mean judgments across given elasticity - $target_id"
        filename = "average_judgment_per_elasticity_"

        sim_mean_elasticity = @chain sim_data begin
            @groupby :gtElasticity
            @combine :judgment = mean(:judgment)
        end

        human_mean_elasticity = @chain human_data begin
            @groupby :gtElasticity
            @combine :judgment = mean(:judgment)
        end
    end

    # @infiltrate
    p = palette(:jet)
    scatter(sim_data.judgment,
            human_data.judgment,
            yerror=human_data.std_err_mean ./ 2,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor=human_mean_elasticity.gtElasticity,
            xlims=(0, 1.05),
            ylims=(0, 1),
            clims=(0, 1),
            xlabel="Model",
            ylabel="Human",
            title=title,
            legend=false,
            colorbar=true,
            palette=p,
            markershape=marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(sim_mean_elasticity.judgment, human_mean_elasticity.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)

    plots_path = generate_plot_path(expt_id, model_id, target_id, noise_id)
    savefig(joinpath(plots_path, string(filename, target_id, "Exp", expt_id, ".png")))
    # gui()

end


function main()

    args = Args()

    # Generate output filepath
    noise_id = generate_noise_id(args)
    plots_path = generate_plot_path(args.expt_id, args.model_id, args.target_id, noise_id)
    println("PLOTS PATH: ", plots_path)

    # Read simulation data and organize into desired form
    sim_data = read_simulation_data(args.expt_id, args.model_id, args.target_id, noise_id, args.algorithm)
    sim_data, error_trials = process_individual_stimuli_sim(sim_data)
    display(sim_data)
    println("High Error Trials: ", error_trials)

    # Plotting variables
    if args.target_id == "Cube"
        marker_shape = :square
    else
        marker_shape = :circle
    end

    # Plot inferred elasticity against ground truth, display correlation
    plot_sim_vs_gt(sim_data, marker_shape, plots_path, 1, args.target_id)

    #=
    # Read human data
    human_data = read_subject_data(1)
    human_data = process_individual_stimuli_human(args.expt_id)

    for expt_id = 1:2
        plot_sim_vs_human(expt_id, target_id, false)
        plot_sim_vs_human(expt_id, target_id, true)
    end
    =#

end

main()
