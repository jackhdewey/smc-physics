# Data analysis
# Plots model-inferred elasticities against human judgments
#
# TODO: Compute error to identify trials with poor model / human correlation

using Plots
using DataFrames
using DataFramesMeta
using Statistics
using Gen
using CSV
using Plots

project_path = dirname(@__DIR__)
println(project_path)

model_id = "Modelv5"
target_id = "Sphere"
var_id = "PosVar075"

# for target_id in ["Cube", "Sphere"]
for expt = 1:2
    plot_individual_stimuli_judgments(expt, target_id)
end

#=
for target_id in ["Cube", "Sphere"]
    for expt = 1:4
        plot_human_vs_gt(expt)
        plot_sim_vs_gt(expt, target_id)
        if expt <= 2
            plot_mean_elasticity_judgments(expt, target_id)
        end
    end
end
=#


# Read the entire folder of simulation data
function read_simulation_data(expt, target_id)
    
    simulation_folder = joinpath(project_path, "Data", "BulletData", model_id, target_id, var_id, "Exp" * string(expt), "Inferences")
    for file in readdir(simulation_folder)
        full_file_path = joinpath(simulation_folder, file)
        data = CSV.read(full_file_path, DataFrame)

        # unpack filename
        _, elasticity_string, variation = split(fname, ['_', '.'])
        elasticity = parse(Int, elasticity_string[end]) * 0.1 # fix order of magnitude
    
        data = insertcols(
            data,
            "filename" => fname,
            "stimulusID" => target_id * "_" * elasticity_string * "_" * variation,
            "gtElasticity" => elasticity,
            "variation" => parse(Int64, variation[4:end])
        )
        push!(all_data, data)
    end

    all_data = vcat(all_data...)    # join all dfs
    # filter(:variation => x -> x <= 15, all_data)
    
    return all_data
end

# For each stimulus, 
function process_individual_stimuli_sim(expt, target_id)
    sim_data = read_simulation_data(expt, target_id)
    sim_data_pred = @chain sim_data begin
        @groupby :stimulusID
        @combine begin
            :judgment = mean(:elasticity)       # :elasticity = :gtElasticity
            :elasticity = first(:gtElasticity)
        end
        # @subset :elasticity .> 0.6
        @orderby :stimulusID
    end

    return sim_data_pred
end

#
function read_gt_data(expt)
    sub_data = read_subject_data(expt)
    gt_data = @select(sub_data, :trialID, :elasticity, :trialType, :stimulusID)
    return gt_data
end

# Read the human predictions
function read_subject_data(expt)

    folders = Dict(
        1 => "Exp1_allElasticities_fullMotion",
        2 => "Exp2_allElasticities_1second",
        3 => "Exp3_mediumElasticity_fullMotion",
        4 => "Exp4_mediumElasticity_1second"
    )
    exp_data_folder = joinpath(project_path, "Data", "HumanData", "EstimationTask", folders[expt], "Results")

    data = []
    for fname in readdir(exp_data_folder)
        sub_data = CSV.read(joinpath(exp_data_folder, fname), DataFrame)
        sub_data = insertcols(sub_data, :filename => fname)
        push!(data, sub_data) # this is a vector of dataframes
    end
    data_df = vcat(data...)     # combine all together with splat operator

    # elasticity is coded as integer
    if expt >= 3
        data_df.elasticity = data_df.elasticity / 10
    end

    return data_df
end

#
function process_individual_stimuli_human(expt)
    sub_data = read_subject_data(expt)
    nsubs = length(unique(sub_data.filename))
    sub_data_pred = @chain sub_data begin
        @groupby :stimulusID
        # @groupby :elasticity# :filename
        # @DataFramesMeta.transform :prediction = mean(:rating) # Gen also has a transform macro
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


#########
# PLOTS #
#########

function generate_plot_path(expt)
    plots_path = joinpath(project_path, "Analysis", "Plots", model_id, target_id, var_id, "Exp" * string(expt), "Judgments")
    if !isdir(plots_path)
        mkdir(plots_path)
    end
    return plots_path
end

#
function plot_sim_vs_gt(expt, target_id)
    sim_raw = read_simulation_data(expt, target_id)

    sim = @chain sim_raw begin
        @groupby :stimulusID
        @combine begin
            :estimate = mean(:elasticity) # :elasticity = :gtElasticity
            :gtElasticity = first(:gtElasticity)
            :std_err_mean = std(:elasticity)
        end
    end

    p = palette(:jet)
    # default(aspect_ratio = :equal)
    scatter(sim.gtElasticity,
            sim.estimate,
            yerror=sim.std_err_mean ./ 2,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor=sim.gtElasticity,# zcolor = :gtElasticity,
            clims=(0, 1),
            xlims=(0, 1.05),
            ylims=(0, 1),
            xlabel="Ground Truth Elasticity",
            ylabel="Model Estimate",
            title="Exp $expt stimulus-specific elasticity ratings - $sim_object",
            colorbar=true,
            legend=false,
            palette=p,
            markershape=marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(sim.gtElasticity, sim.estimate), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)
    plots_path = generate_plot_path(expt)
    savefig(joinpath(plots_path, string("individual_stimuli_judgments_", "against_gt_model_", target_id, "Exp", expt, ".png")))
end

#
function plot_human_vs_gt(expt)
    human = process_individual_stimuli_human(expt)
    # scatter(human.gtElasticity, human.judgment)
    # @autoinfiltrate
    p = palette(:jet)
    # default(aspect_ratio = :equal)
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
    plots_path = generate_plot_path(expt)
    savefig(joinpath(plots_path, string("individual_stimuli_judgments_", "against_gt_human", "Exp", expt, ".png")))

    # gui()
end

#
function plot_individual_stimuli_judgments(expt, target_id)
    human = process_individual_stimuli_human(expt)
    sim = process_individual_stimuli_sim(expt, target_id)

    p = palette(:jet)
    # default(aspect_ratio = :equal)
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
    plots_path = generate_plot_path(expt)
    savefig(joinpath(plots_path, string("individual_stimuli_judgments_high_elasticity", target_id, "Exp", expt, ".png")))
    gui()
end

#
function plot_mean_elasticity_judgments(expt, target_id)

    human = process_individual_stimuli_human(expt)
    sim = process_individual_stimuli_sim(expt, target_id)
    p = palette(:jet)
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
    plots_path = generate_plot_path(expt)
    savefig(joinpath(plots_path, string("average_judgment_per_elasticity_", target_id, "Exp", expt, ".png")))
    # gui()
end
