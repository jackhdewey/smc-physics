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

sim_object = "Cube"

project_path = dirname(@__DIR__)
folder = joinpath(project_path, "BulletData", sim_object, "Observations/")

function read_simulation_file(fname)
    data = CSV.read(fname, DataFrame)

    # unpack filename
    sim_object, elasticity_string, replicate = split(fname, ['_', '.'])
    elasticity = parse(Int, elasticity_string[end]) * 0.1 # fix order of magnitude
    return data, elasticity
end

function get_particle_elasticity_estimate(data, i)
    # get particle i
    particle_vals = filter(:particle => x -> x == i, data)

    # because the elasticity values are all the same for a given particle in a given file
    row = first(particle_vals)
    return row.elasticity
end

function calc_average_elasticity_estimate()
    elasticity_estimates = Dict()

    for file in readdir(folder)
        full_file_path = joinpath(folder, file)

        data, elasticity = read_simulation_file(full_file_path)
        estimates_this_data_file = []

        for i = 1:maximum(data.particle) # total # particles
            # get estimate for particle i
            particle_estimate = get_particle_elasticity_estimate(data, i)

            push!(estimates_this_data_file, particle_estimate)
        end

        if elasticity in keys(elasticity_estimates)
            push!(elasticity_estimates[elasticity], mean(estimates_this_data_file))
        else
            elasticity_estimates[elasticity] = [mean(estimates_this_data_file)]
        end
    end
    mean_particle_elasticity = Dict((i) => mean(elasticity_estimates[i]) for i in keys(elasticity_estimates))
    mean_particle_elasticity = sort(mean_particle_elasticity)

    errorbars = Dict((i) => std(elasticity_estimates[i]) / length(elasticity_estimates[i]) for i in keys(elasticity_estimates))
    errorbars = sort(errorbars)

    return (mean_particle_elasticity, errorbars)
end

function plot_all_estimates()
    all_estimates = []
    for file in readdir(folder)
        full_file_path = joinpath(folder, file)
        data, rf_ela = read_simulation_file(full_file_path)

        for i = 1:maximum(data.particle)
            est_ela = get_particle_elasticity_estimate(data, i)
            push!(all_estimates, [rf_ela, est_ela])
        end
    end

    ests_mat = hcat(all_estimates...)'

    scatter(
        ests_mat[:, 1],
        ests_mat[:, 2],
        markeralpha=0.1,
        markersize=0.2,
        xlims=(0, 1),
        ylims=(0, 1),
        aspect_ratio=:equal,
        xlabel="Ground truth elasticity",
        ylabel="Estimated elasticity",
        title="All " * sim_object * " Simulations",
        legend=false
    )
end

function plot_individual_stimuli()

        full_file_path = vcat([CSV.read(joinpath(folder, file), DataFrame) for file in readdir(folder)])

        data, elasticity = read_simulation_file(full_file_path)

end

function calc_model_predictions()


    include(joinpath(project_path, "Utilities/fileio.jl"))

    mean_elasticity, errorbars = calc_average_elasticity_estimate()

    rf = collect(keys(mean_elasticity))
    estimates = collect(values(mean_elasticity))
    errorbars = collect(values(errorbars))

    plot(
        rf,
        estimates,
        yerror=errorbars,
        xlabel="Ground truth elasticity",
        ylabel="Estimated elasticity",
        xlims=(0, 1),
        ylims=(0, 1),
        aspect_ratio=:equal,
        title="Averaged Estimates for " * sim_object,
        legend=false
    )

    corr_string = "r = " * string(round(cor(estimates, rf), digits=3))
    annotate!(0.5, 0.3, corr_string)

    savefig("average_elasticity_estimates_" * sim_object * ".png")

    plot_all_estimates()
    savefig("all_elasticity_estimates_" * sim_object * ".png")

    # println("Correlation between estimates and ground truth for " * sim_object * ": " * string(cor(estimates, rf)))
    estimates
end

function read_subject_data(path)
    data = []
    for fname in readdir(path)
        println(fname)
        sub_data = CSV.read(joinpath(project_path, path, fname), DataFrame)
        sub_data = insertcols(sub_data, :ID => fname)
        data = vcat(data, sub_data) # this is a vector of dataframes
    end
    data = vcat(data...)        # combine all together with splat operator

end

function get_sub_predictions()
    sub_data = read_subject_data(joinpath(project_path, "SubjectDataCubes"))

    preds = @chain sub_data begin
        @groupby(:elasticity)
        @combine(:prediction = mean(:rating))
        @orderby(:elasticity)
    end
    preds.prediction
end

sub_predictions = get_sub_predictions()
model_predictions = calc_model_predictions()

scatter(model_predictions, sub_predictions, aspect_ratio=:equal, legend=false, xlabel="Model elasticity",
ylabel = "Subject prediction", title = sim_object)
plot!(0:1, 0:1)
gui()
savefig(string("sub_model_", sim_object, ".png"))
cor(model_predictions, sub_predictions)

# plot_individual_stimuli()
