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

# Group subject elasticity judgments by stimulus and compute mean rating for each stimulus
function process_individual_stimuli_human(expt_id)

    sub_data = read_subject_data(expt_id)
    nsubs = length(unique(sub_data.filename))

    sub_data_pred = @chain sub_data begin
        @groupby :stimulusID
        # @groupby :elasticity# :filename
        # @DataFramesMeta.transform :prediction = mean(:rating)     # Gen also has a transform macro
        @combine begin
            :judgment = mean(:rating)
            :std_err_mean = std(:rating) / sqrt(nsubs) # number of subjects
            :gtElasticity = mean(:elasticity)
        end
        # @subset :gtElasticity .> 0.6
        @orderby :stimulusID
    end

    return sub_data_pred
end

# Compute error on each stimulus filter to extract high error trials
function process_individual_stimuli_sim(expt_id, model_id, target_id, noise_id, alg)

    # Define model's elasticity judgment as average of output particles
    sim_data = read_simulation_data(expt_id, model_id, target_id, noise_id, alg)
    display(sim_data)
    sim_data_pred = @chain sim_data begin
        @groupby :stimulusID
        @combine begin

            :judgment = mean(:elasticity)       
            :elasticity = first(:gtElasticity)  # :elasticity = :gtElasticity

            # compute error w.r.t. ground truth
            :error = mean(:elasticity) - first(:gtElasticity)
        end
        # @subset :elasticity .> 0.6
        @orderby :stimulusID
    end

    # Filter by error threshold, save high error trials
    high_error = sim_data_pred[sim_data_pred.error .> .2, :]
    high_error_trials = high_error[:, 1]
    println(high_error_trials)

    return sim_data_pred, high_error_trials
end


#########
# PLOTS #
#########

# Plot model against ground truth
function plot_sim_vs_gt(sim, expt_id, target_id, marker_shape, plots_path)

    # default(aspect_ratio = :equal)
    p = palette(:jet)
    scatter(sim.elasticity,
            sim.judgment,
            #yerror=sim.std_err_mean ./ 2,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor=sim.elasticity,          # zcolor = :gtElasticity,
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
    corr_string = "r = " * string(round(cor(sim.elasticity, sim.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)
    savefig(joinpath(plots_path, string("individual_stimuli_judgments_", "against_gt_model_", target_id, "Exp", expt_id, ".png")))

end

# Plot human judgments against the ground truth
function plot_human_vs_gt(expt_id)
    
    human = process_individual_stimuli_human(expt_id)

    # default(aspect_ratio = :equal)
    p = palette(:jet)
    scatter(human.gtElasticity,
        human.judgment,
        yerror=human.std_err_mean ./ 2,
        aspect_ratio=:equal,
        markersize=5,
        markeralpha=0.5,
        zcolor=human.gtElasticity,# zcolor = :gtElasticity,
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
    corr_string = "r = " * string(round(cor(human.gtElasticity, human.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)
    plots_path = generate_plot_path(expt_id, model_id, target_id, noise_id)
    savefig(joinpath(plots_path, string("individual_stimuli_judgments_", "against_gt_human", "Exp", expt_id, ".png")))
    # gui()

end


# Plot mean model estimates against mean human judgments
function plot_sim_vs_human(expt_id, target_id, plot_mean)

    human = process_individual_stimuli_human(expt_id)
    sim = process_individual_stimuli_sim(expt_id, model_id, target_id, noise_id)

    # default(aspect_ratio = :equal)

    title = "Exp $expt_id stimulus-specific elasticity ratings - $target_id"
    filename = "individual_stimuli_judgments_high_elasticity"

    if plot_mean

        title = "Exp $expt_id mean judgments across given elasticity - $target_id"
        filename = "average_judgment_per_elasticity_"

        sim = @chain sim begin
            @groupby :elasticity
            @combine :judgment = mean(:judgment)
        end

        human = @chain human begin
            @groupby :gtElasticity
            @combine :judgment = mean(:judgment)
        end
    end

    # @infiltrate
    p = palette(:jet)
    scatter(sim.judgment,
            human.judgment,
            yerror=human.std_err_mean ./ 2,
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
    corr_string = "r = " * string(round(cor(model_mean_elasticity.judgment, human_mean_elasticity.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)
    plots_path = generate_plot_path(expt_id, model_id, target_id, noise_id)
    savefig(joinpath(plots_path, string(filename, target_id, "Exp", expt_id, ".png")))
    # gui()

end


function main()

    args = Args()

    # Generate output filepath
    noise_id = generate_noise_id(args)
    sim, _ = process_individual_stimuli_sim(args.expt_id, args.model_id, args.target_id, noise_id, args.algorithm)
    plots_path = generate_plot_path(args.expt_id, args.model_id, args.target_id, noise_id)
    println("PLOTS PATH: ", plots_path)

    # Plotting variables
    if args.target_id == "Cube"
        marker_shape = :square
    else
        marker_shape = :circle
    end

    #process_individual_stimuli_sim(args.expt_id, args.model_id, args.target_id, args.noise_id)
    #plot_sim_vs_human_individual_stimuli(model_id, target_id, noise_id, expt_id)
    plot_sim_vs_gt(sim, 1, args.target_id, marker_shape, plots_path)

    #=
    for expt_id = 1:2
        plot_sim_vs_human(expt_id, target_id, false)
        plot_sim_vs_human(expt_id, target_id, true)
    end
    =#

end

main()
