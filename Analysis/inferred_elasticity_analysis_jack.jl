# Data analysis
# Assess model performance at inferring elasticity, relative to ground truth and/or human judgments
# Identify trials with poor model / human correlation by absolute error
# Generate 2D plots showing model-inferred elasticities vs human judgments

using DataFrames
using DataFramesMeta
using Statistics
using Plots

include("utils.jl")

project_path = dirname(@__DIR__)
println(project_path)

model_id = "Modelv5"
target_id = "Cube"
noise_id = "PosVar05"
expt_id = "BulletxBullet"

if target_id == "Cube"
    marker_shape = :square
else
    marker_shape = :circle
end

#################
# Data Analysis #
#################

# Compute particle filter error on each stimulus, filter to extract high error trials
function process_individual_stimuli_sim(model_id, target_id, noise_id, expt_id)

    # Define model's elasticity judgment as average of output particles, compute error w.r.t. ground truth
    sim_data = read_simulation_data(model_id, target_id, noise_id, expt_id)
    sim_data_pred = @chain sim_data begin
        @groupby :stimulusID
        @combine begin
            :judgment = mean(:elasticity)       # :elasticity = :gtElasticity
            :elasticity = first(:gtElasticity)
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

# Group subject elasticity judgments by stimulus and compute mean rating for each stimulus
function process_individual_stimuli_human(expt)

    sub_data = read_subject_data(expt)
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

    println(sub_data_pred)

    return sub_data_pred
end


#########
# PLOTS #
#########

# Plot model against ground truth
function plot_sim_vs_gt(expt, target_id)

    sim, _ = process_individual_stimuli_sim(model_id, target_id, noise_id, expt_id)

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
            title="Exp $expt stimulus-specific elasticity ratings - $target_id",
            colorbar=true,
            legend=false,
            palette=p,
            markershape=marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(sim.elasticity, sim.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)
    plots_path = generate_plot_path(model_id, target_id, noise_id, expt_id)
    savefig(joinpath(plots_path, string("individual_stimuli_judgments_", "against_gt_model_", target_id, "Exp", expt, ".png")))

end

# Plot human judgments against the ground truth
function plot_human_vs_gt(expt)
    
    human = process_individual_stimuli_human(expt)

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
        title="Exp $expt stimulus-specific human ratings",
        colorbar=true,
        legend=false,
        palette=p,
        markershape=marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(human.gtElasticity, human.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)
    plots_path = generate_plot_path(model_id, target_id, noise_id, expt_id)
    savefig(joinpath(plots_path, string("individual_stimuli_judgments_", "against_gt_human", "Exp", expt, ".png")))
    # gui()

end

# Plot model estimates against human judgments
function plot_sim_vs_human_individual_stimuli(expt, target_id)

    human = process_individual_stimuli_human(expt)
    sim = process_individual_stimuli_sim(model_id, target_id, noise_id, expt_id)

    # default(aspect_ratio = :equal)
    p = palette(:jet)
    scatter(sim.judgment,
            human.judgment,
            yerror=human.std_err_mean ./ 2,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor=human.gtElasticity,# zcolor = :gtElasticity,
            clims=(0, 1),
            xlims=(0, 1.05),
            ylims=(0, 1),
            xlabel="Model",
            ylabel="Human",
            title="Exp $expt stimulus-specific elasticity ratings - $target_id",
            colorbar=true,
            legend=false,
            palette=p,
            markershape=marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(sim.judgment, human.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)
    plots_path = generate_plot_path(model_id, target_id, noise_id, expt_id)
    savefig(joinpath(plots_path, string("individual_stimuli_judgments_high_elasticity", target_id, "Exp", expt, ".png")))
    gui()

end

# Plot mean model estimates against mean human judgments
function plot_sim_vs_human_mean(expt, target_id)

    human = process_individual_stimuli_human(expt)
    sim = process_individual_stimuli_sim(model_id, target_id, noise_id, expt_id)

    # default(aspect_ratio = :equal)
    model_mean_elasticity = @chain sim begin
        @groupby :elasticity
        @combine :judgment = mean(:judgment)
    end

    human_mean_elasticity = @chain human begin
        @groupby :gtElasticity
        @combine :judgment = mean(:judgment)
    end

    # @infiltrate
    p = palette(:jet)
    scatter(model_mean_elasticity.judgment,
            human_mean_elasticity.judgment,
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
            title="Exp $expt mean judgments across given elasticity - $target_id",
            # legend=false,
            legend=false,
            colorbar=true,
            palette=p,
            markershape=marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(model_mean_elasticity.judgment, human_mean_elasticity.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)
    plots_path = generate_plot_path(model_id, target_id, noise_id, expt_id)
    savefig(joinpath(plots_path, string("average_judgment_per_elasticity_", target_id, "Exp", expt, ".png")))
    # gui()

end

#process_individual_stimuli_sim(model_id, target_id, noise_id, expt_id)
#plot_sim_vs_human_individual_stimuli(model_id, target_id, noise_id, expt_id)
plot_sim_vs_gt(1, target_id)

#=
for expt = 1:2
    plot_sim_vs_human_individual_stimuli(expt, target_id)
    plot_sim_vs_human_mean(expt, target_id)
end
=#
