# TODO: Visualize by plotting bounce locations in 3D

include("Utilities/plots.jl")
include("Utilities/fileio.jl")


function main()

    all_files = filter_unwanted_filenames(readdir("BulletData/Intermediate/")) 
    sort!(all_files, lt=trial_particle_order)

    for file in all_files
        tokens = split(file, "_")
        current_trial = tokens[3]
        for i=1:30
            particles = CSV.read(string("BulletData/Intermediate/", file), DataFrame)
            plot_traj(particles)
        end
    end
end

main()

