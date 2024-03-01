include("Utilities/plots.jl")

function main()

    all_files = readdir("BulletData/Intermediate/")

    for file in all_files
        tokens = split(file, "_")
        print(file, "\n")
        #all_particles = CSV.read(file, DataFrame)
        #plot_traj(all_particles)
    end
end

main()

