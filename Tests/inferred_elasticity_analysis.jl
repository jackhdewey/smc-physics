using DataFrames
using Statistics
using Gen
using CSV
using Plots

# TODO more ground truth RF sphere data
# TODO different elasticity of cube data (later?)
# TODO look at trajectories of sphere data
# TODO look at stats of sphere data

project_path = dirname(@__DIR__)
folder = joinpath(project_path, "BulletData/Sphere/Observations/")

include(joinpath(project_path, "Utilities/fileio.jl"))


function read_file(fname)
    data = CSV.read(fname, DataFrame)

    # unpack filename
    object, elasticity_string, replicate = split(fname, ['_', '.'])
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

        data, elasticity = read_file(full_file_path)
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
        data, rf_ela = read_file(full_file_path)

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
        xlabel="Ground truth elasticity",
        ylabel="Estimated elasticity",
        title="All simulations",
        legend=false
    )
end

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
    title="Averages",
    legend=false
)
savefig("average.png")

plot_all_estimates()
savefig("all_estimates.png")
