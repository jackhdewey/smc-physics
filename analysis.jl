# Data analysis module

# Organization: grid by trial, with 30 grids for each time step, gris is 1x1x1
# TODO: Set the alpha / intensity to reflect the log weight
# TODO: Add the ground truth to each plot

include("Utilities/plots.jl")
include("Utilities/fileio.jl")


function main()

    plt1 = plot3d( 
        1,
        xlim = (-.5, .5),
        ylim = (-.5, .5),
        zlim = (0, 1),
        title = "Particle Filter Flow",
        legend = false,
        marker = 2,
        seriestype=:scatter
    )
    plt2 = plot3d( 
        1,
        xlim = (-.5, .5),
        ylim = (-.5, .5),
        zlim = (0, 1),
        title = "Particle Filter Flow",
        legend = false,
        marker = 2,
        seriestype=:scatter
    )

    all_files = filter_unwanted_filenames(readdir("BulletData/Intermediate/")) 

    sort!(all_files, lt=trial_particle_order)
    file = all_files[30]

    #for file in all_files

        data = CSV.read(string("BulletData/Intermediate/", file), DataFrame)

        ground_truth = CSV.read("RealFlowData/Cube_Ela3_Var26_observed.csv", DataFrame)

        true_x = []
        true_y = []
        true_z = []
        for i=1:30
            true_x = [true_x; ground_truth[i, 1]]
            true_z = [true_z; ground_truth[i, 2]]
            true_y = [true_y; ground_truth[i, 3]]
        end

        for i=1:20
            display(data[data.particle .== i, 5:7])
            elasticity = data[data.particle .== i, 2]
            weight = data[data.particle .== i, 3]
            print(string("Elasticity: ", elasticity[1], "\n"))
            print(string("Log Weight: ", weight[1], "\n"))
        end

        tokens = split(file, "_")
        time_step = parse(Int64, replace(tokens[4], ".csv" => ""))
        
        # Iterate over each particle
        for i=1:time_step:20*time_step            

            x_trajectory = []
            y_trajectory = []
            z_trajectory = []

            # Iterate over each time step
            for j=1:time_step
                row = (1) + j
                x_trajectory = [x_trajectory; data[row, 5]]
                y_trajectory = [y_trajectory; data[row, 6]]
                z_trajectory = [z_trajectory; data[row, 7]]
            end

            # Plot the 
            for k=1:length(x_trajectory)
                push!(plt1, x_trajectory[k], y_trajectory[k], z_trajectory[k])
            end

            for k=1:length(true_x)
                push!(plt2, true_x[k], true_y[k], true_z[k])
            end
        end
        display(plt1)
        display(plt2)
    #end
end

main()

