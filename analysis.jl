# Data analysis module

# Organization: grid by trial, with 30 grids for each time step, gris is 1x1x1
# TODO: Set the alpha / intensity to reflect the log weight
# TODO: Add the ground truth to each plot


include("Utilities/plots.jl")
include("Utilities/fileio.jl")


function main()

    all_files = filter_unwanted_filenames(readdir("BulletData/Intermediate/")) 
    sort!(all_files, lt=trial_particle_order)
    file = all_files[2]

    #for file in all_files

        data = CSV.read(string("BulletData/Intermediate/", file), DataFrame)

        tokens = split(file, "_")
        time_step = parse(Int64, replace(tokens[4], ".csv" => ""))
        
        trajectories = cat([], dims = 3)
        for i=1:time_step:20*time_step
            trajectory = []
            for j=1:time_step
                col = (i-1) + j
                x_pos = data[col, 5]
                y_pos = data[col, 6]
                z_pos = data[col, 7]
                trajectory = [trajectory; [x_pos, y_pos, z_pos]]
            end
            trajectories = [trajectories; trajectory]
        end

        println(trajectories)
    #end
    #plt = plot3d()
end

main()

