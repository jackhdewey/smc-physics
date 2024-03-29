# Data analysis module
# Organization: 30 plots for each triel, one for each time step - ground truth trajectory plotted in red

# TODO: Set the alpha / intensity to reflect the log weight
# TODO: Save plots in a folder
# TODO: Allow more flexible selection of elasticity / trial interval

include("Utilities/plots.jl")
include("Utilities/fileio.jl")


function main()

    for i=1:10

        # Generate plot base
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

        trial = i

        # Read ground truth trajectory
        gt_file = string("RealFlowData/Sphere/SoftSphere_Ela3_Var", trial, "_observed.csv")
        print(gt_file)
        ground_truth = CSV.read(gt_file, DataFrame)
        true_x = []
        true_y = []
        true_z = []

        # Read all particle files
        all_files = filter_unwanted_filenames(readdir("BulletData/Sphere/Intermediate/")) 
        sort!(all_files, lt=trial_particle_order)

        # For each time step in the first trial
        for x=1:30

            index = 30 * (trial-1) + x

            true_x = [true_x; ground_truth[x, 1]]
            true_y = [true_y; ground_truth[x, 2]]
            true_z = [true_z; ground_truth[x, 3]]

            plot!(plt1, true_x, true_y, true_z, linewidth=3, linecolor=:red)

            file = all_files[index]
            data = CSV.read(string("BulletData/Sphere/Intermediate/", file), DataFrame)
            tokens = split(file, "_")
            time_step = parse(Int64, replace(tokens[4], ".csv" => ""))
        
            # For each particle at this time step
            for i=1:20    

                particle = data[data.particle .== i, :]

                elasticity = data[data.particle .== i, 2]
                weight = data[data.particle .== i, 3]

                # @df plot!(plt1, particle[:, 5:7])

                # Iterate over each time step and push to plot
                x_trajectory = []
                y_trajectory = []
                z_trajectory = []
                for row=1:time_step
                    x_trajectory = [x_trajectory; particle[row, 5]]
                    y_trajectory = [y_trajectory; particle[row, 6]]
                    z_trajectory = [z_trajectory; particle[row, 7]]
                end

                plot!(plt1, x_trajectory, y_trajectory, z_trajectory)

            end

            # Plot the current particles
            display(plt1)
        end
    end
end

main()
