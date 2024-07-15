using Plots
using DataFrames
using DataFramesMeta
using Statistics
using Gen
using CSV
using Plots

project_path = dirname(@__DIR__)



# for variable in ["ObsVar025", "ObsVar075", "PosVar025", "PosVar05", "PosVar075"]
for variable in ["PosVar075"]

    sim_object = "Sphere"

    if sim_object == "Cube"
        marker_shape = :square
    else
        marker_shape = :circle
    end

    plots_path = joinpath(project_path, "Analysis", "Plots", "Modelv5", variable, "Exp1", sim_object)
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    function read_simulation_file(fname, sim_object)
        data = CSV.read(fname, DataFrame)

        # unpack filename
        _, elasticity_string, variation = split(fname, ['_', '.'])
        elasticity = parse(Int, elasticity_string[end]) * 0.1 # fix order of magnitude

        data = insertcols(
            data,
            "filename" => fname,
            "stimulusID" => sim_object * "_" * elasticity_string * "_" * variation,
            "gtElasticity" => elasticity,
            "variation" => parse(Int64, variation[4:end])
        )
        return data
    end

    function read_simulation_data(expt, sim_object)
        # simulation_folder = joinpath(project_path, "BulletData", "Modelv5", "Exp" * string(expt), sim_object, "Inferences")
        simulation_folder = joinpath(project_path, "Data", "BulletData", "Modelv5", "Sphere", variable, "Exp1", "Inferences")
        all_data = []
        for file in readdir(simulation_folder)
            full_file_path = joinpath(simulation_folder, file)
            data = read_simulation_file(full_file_path, sim_object)
            push!(all_data, data)
        end
        all_data = vcat(all_data...)    # join all dfs
        # return filter(:variation => x -> x <= 15, all_data)
        return all_data
    end

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

    function read_gt_data(expt)

        sub_data = read_subject_data(expt)
        gt_data = @select(sub_data, :trialID, :elasticity, :trialType, :stimulusID)
        return gt_data
    end

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
        savefig(joinpath(plots_path, string("individual_stimuli_judgments_", "against_gt_human", "Exp", expt, "_", variable, ".png")))

        # gui()
    end

    function plot_sim_vs_gt(expt, sim_object)
        sim_raw = read_simulation_data(expt, sim_object)

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
        savefig(joinpath(plots_path, string("individual_stimuli_judgments_", "against_gt_model_", sim_object, "Exp", expt, ".png")))
    end

    function process_individual_stimuli_sim(expt, sim_object)
        sim_data = read_simulation_data(expt, sim_object)
        sim_data_pred = @chain sim_data begin
            @groupby :stimulusID
            @combine begin
                :judgment = mean(:elasticity) # :elasticity = :gtElasticity
                :elasticity = first(:gtElasticity)
            end
            # @subset :elasticity .> 0.6
            @orderby :stimulusID
        end

        print(size(sim_data_pred))
        return sim_data_pred
    end

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

    function plot_individual_stimuli_judgments(expt, sim_object)
        human = process_individual_stimuli_human(expt)
        sim = process_individual_stimuli_sim(expt, sim_object)

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
            title="Exp $expt stimulus-specific ratings - $sim_object - $variable",
            colorbar=true,
            legend=false,
            palette=p,
            markershape=marker_shape)

        plot!(0:1, 0:1, line=:dash)
        corr_string = "r = " * string(round(cor(sim.judgment, human.judgment), digits=3))
        annotate!(0.2, 0.8, corr_string, 10)
        savefig(joinpath(plots_path, string("individual_stimuli_judgments_high_elasticity", sim_object, "Exp", expt, variable, ".png")))
        gui()
    end

    function plot_mean_elasticity_judgments(expt, sim_object)

        human = process_individual_stimuli_human(expt)
        sim = process_individual_stimuli_sim(expt, sim_object)
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
            title="Exp $expt mean judgments for elasticity - $sim_object - $variable",
            # legend=false,
            legend=false,
            colorbar=true,
            palette=p,
            markershape=marker_shape)


        plot!(0:1, 0:1, line=:dash)
        corr_string = "r = " * string(round(cor(model_mean_elasticity.judgment, human_mean_elasticity.judgment), digits=3))
        annotate!(0.2, 0.8, corr_string, 10)
        savefig(joinpath(plots_path, string("average_judgment_per_elasticity_", sim_object, "Exp", expt, "_", variable, ".png")))
        # gui()
    end

    # for sim_object in ["Cube", "Sphere"]
    for expt = 1:1
        plot_individual_stimuli_judgments(expt, sim_object)
        # if expt <= 2
        plot_mean_elasticity_judgments(expt, sim_object)
        # end
    end
    # end
end

# for sim_object in ["Cube"# , "Sphere"
#                    ]
#     for expt = 1:4
#         plot_human_vs_gt(expt)
#         plot_sim_vs_gt(expt, sim_object)
#         # if expt <= 2
#         #     plot_mean_elasticity_judgments(expt, sim_object)
#         # end
#     end
# end
