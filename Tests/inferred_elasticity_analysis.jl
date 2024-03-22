using DataFrames
using Statistics
using Gen
using CSV
# using Plots

# TODO more ground truth RF sphere data
# TODO different elasticity of cube data (later?)
# TODO look at trajectories of sphere data
# TODO look at stats of sphere data

project_path = dirname(@__DIR__)
include(joinpath(project_path, "Utilities/fileio.jl"))

function main()

    elasticity_estimates = Dict()

    folder = joinpath(project_path, "BulletData/Cube/Observations/")
    for file in readdir(folder)
        full_file_path = joinpath(folder, file)

        # obs = read_intermediate_file(full_file_path)
        data = CSV.read(full_file_path, DataFrame)
        # unpack filename
        _, elasticity, replicate = split(file, ['_', '.'])

        estimates_this_data_file = []
        for i = 1:maximum(data.particle) # total # particles
            # get particle i
            particle_vals = filter(:particle => x -> x == i, data)
            # because the elasticity values are all the same for a given particle in a given file
            row = first(particle_vals)

            push!(estimates_this_data_file, row.elasticity)

        end

        if elasticity in keys(elasticity_estimates)
            push!(elasticity_estimates[elasticity], mean(estimates_this_data_file))
        else
            elasticity_estimates[elasticity] = [mean(estimates_this_data_file)]
        end
    end
    mean_particle_elasticity = Dict(i => mean(elasticity_estimates[i]) for i in keys(elasticity_estimates))
end

data = sort(main())
