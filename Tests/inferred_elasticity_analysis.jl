using Plots: StatsBase, get_aspect_ratio
using DataFrames
using DataFramesMeta
using Statistics
using Gen
using CSV
using Plots

# TODO more ground truth RF sphere data
# TODO different elasticity of cube data (later?)
# TODO look at trajectories of sphere data
# TODO look at stats of sphere data


# plot defaults

sim_object = "Cube"

project_path = dirname(@__DIR__)
simulation_folder = joinpath(project_path, "BulletData", sim_object, "Variable Frames", "Exp1", "Observations")
# simulation_folder = joinpath(project_path, "BulletData", sim_object, "30 Frames", "Observations")

subject_folder = "/Users/maxs/smc-physics/Data&StimuliForGenModel/"

function read_simulation_file(fname)
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

function read_all_simulation_data()
    all_data = []
    for file in readdir(simulation_folder)
        full_file_path = joinpath(simulation_folder, file)
        data = read_simulation_file(full_file_path)
        push!(all_data, data)
    end

    all_data = vcat(all_data...)    # join all dfs
    # println(all_data)
    return filter(:variation => x -> x <= 15, all_data)
end



function get_particle_elasticity_estimate(data, i)
    # get particle i
    particle_vals = filter(:particle => x -> x == i, data)

    # because the elasticity values are all the same for a given particle in a given file
    row = first(particle_vals)
    return row.elasticity
end

function get_simulation_elasticity_estimates()
    elasticity_estimates = []

    for file in readdir(simulation_folder)
        full_file_path = joinpath(simulation_folder, file)
        data = read_simulation_file(full_file_path)
        # estimates_this_data_file = []

        # for i = 1:maximum(data.particle) # total # particles
        #     # get estimate for particle i
        #     particle_estimate = get_particle_elasticity_estimate(data, i)
        #     push!(estimates_this_data_file, particle_estimate)
        # vcat!(elasticity_estimates, particle)
        # end

        # @autoinfiltrate
        elasticity_estimates = vcat(data, elasticity_estimates)

    end
    sim_data = vcat(elasticity_estimates...)

    simulation_predictions = @chain sim_data begin
        @groupby
    end

    return
end

function read_subject_data()
    exp_data_folder = joinpath(project_path, "HumanData", "EstimationTask", "Exp1_allElasticities_fullMotion", "Results")
    data = []
    for fname in readdir(exp_data_folder)
        sub_data = CSV.read(joinpath(exp_data_folder, fname), DataFrame)
        sub_data = insertcols(sub_data, :filename => fname)
        # @autoinfiltrate
        # data = vcat(data, sub_data) # this is a vector of dataframes
        # sub_data[!, :rating] = StatsBase.zscore(sub_data[!, :rating])
        push!(data, sub_data) # this is a vector of dataframes
    end

    return vcat(data...)        # combine all together with splat operator
end

function individual_stimuli_sim()
    sim_data = read_all_simulation_data()

    sim_data_pred = @chain sim_data begin
        # @groupby :gtElasticity
        @groupby :stimulusID
        @combine begin
            :judgment = mean(:elasticity) # :elasticity = :gtElasticity
            :elasticity = first(:gtElasticity)
        end
        # @orderby :stimulusID
        # @orderby :gtElasticity
    end

    return sim_data_pred
end
function individual_stimuli_human()
    sub_data = read_subject_data()

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
        @orderby :stimulusID
        # @orderby :elasticity
    end

    # scatter(
    #     sub_data_pred.elasticity,
    #     sub_data_pred.prediction,
    #     aspect_ratio = :equal,
    #     xlims = (0, 1),
    #     ylims = (0, 1),
    #     markeralpha = .2,
    #     # markersize = 2,
    #     xlabel = "Elasticity",
    #     ylabel = "Average Subject Rating",
    #     title = "Stimulus-specific ratings in " * sim_object * " condition",
    #     legend = false
    # )
    # plot!(0:1, 0:1)
    # savefig(string("stimulus_specific_estimates_", sim_object, "_model.png"))
    return sub_data_pred
end

function plot_individual_stimuli_judgments()
    human = individual_stimuli_human()
    sim = individual_stimuli_sim()

    p = Plots.palette(:jet)
    # default(aspect_ratio = :equal)
    scatter(sim.judgment,
            human.judgment,
            yerror=human.std_err_mean ./ 2,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor=human.gtElasticity,# zcolor = :gtElasticity,
            xlims=(0, 1.05),
            ylims=(0, 1),
            xlabel="Model",
            ylabel="Human",
            title="Stimulus-specific Elasticity Ratings",
            legend=false,
            palette=p)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(sim.judgment, human.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)


    gui()
end

function plot_mean_elasticity_judgments()

    human = individual_stimuli_human()
    sim = individual_stimuli_sim()
    print(first(human))
    print(first(sim))
    p = Plots.palette(:jet)
    # default(aspect_ratio = :equal)
    model_mean_elasticity = @chain sim begin
        @groupby :elasticity
        @combine :judgment = mean(:judgment)
    end

    human_mean_elasticity = @chain human begin
        @groupby :gtElasticity
        @combine :judgment = mean(:judgment)
    end

    scatter(model_mean_elasticity.judgment,
            human_mean_elasticity.judgment,
            yerror=human.std_err_mean ./ 2,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor=human.gtElasticity,# zcolor = :gtElasticity,
            xlims=(0, 1.05),
            ylims=(0, 1),
            xlabel="Model",
            ylabel="Human",
            title="Averaged across stimuli",
            legend=false,
            palette=p)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(model_mean_elasticity.judgment, human_mean_elasticity.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)

    gui()
end


# plot_individual_stimuli_judgments()
plot_mean_elasticity_judgments()
# function read_all_subject_data()
#     all_data = []
#     subject_data_folder = joinpath(project_path, "DaHumanData "EstimationTask", "Exp1_allElasticities_fullMotion", "Results")
#     for exp in [dir for dir in readdir(subject_data_folders) if !startswith(dir, ".DS")][1:1]
#         full_file_path = joinpath(subject_data_folders, exp, "Results")
#         println(full_file_path)
#         data = read_subject_data(full_file_path)
#         push!(all_data, data)
#     end

#     return vcat(all_data...)    # join all dfs
# end

# human_judgments = get_sub_judgments()
# model_judgments = calc_model_judgments()

# scatter(
#     model_judgments,
#     human_judgments,
#     aspect_ratio=:equal,
#     legend=false,
#     xlabel="Model elasticity",
#     ylabel="Subject prediction",
#     title=sim_object
# )

# plot!(0:1, 0:1)
# gui()
# savefig(string("sub_model_", sim_object, ".png"))

# cor(model_judgments, human_judgments)

# plot_individual_stimuli()
