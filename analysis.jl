
# TODO: Visualize by plotting bounce locations in 3D

include("Utilities/plots.jl")
include("Utilities/fileio.jl")

function trial_particle_order(x, y)
    x_tokens = split(x, "_")
    y_tokens = split(y, "_")

    trial1 = replace(x_tokens[3], "Var" => "")
    trial2 = replace(y_tokens[3], "Var" => "")

    parse(Int64, trial1)
    parse(Int64, trial2)

    if (trial1 != trial2)
        return trial1 < trial2
    end

    particle1 = replace(x_tokens[4], ".csv" => "")
    particle2 = replace(y_tokens[4], ".csv" => "")

    parse(Int64, particle1)
    parse(Int64, particle2)

    return particle1 < particle2

end

function main()

    all_files = filter_unwanted_filenames(readdir("BulletData/Intermediate/")) 
    sort(all_files, lt=trial_particle_order)

    for file in all_files
        print(file, "\n")
        tokens = split(file, "_")
        #all_particles = CSV.read(file, DataFrame)
        #plot_traj(all_particles)
    end
end

main()

