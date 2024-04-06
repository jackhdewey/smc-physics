# Data analysis module
# Organization: 30 plots for each triel, one for each time step - ground truth trajectory plotted in red

# TODO: Set the alpha / intensity to reflect the log weight
# TODO: Save plots in a folder
# TODO: Allow more flexible selection of elasticity / trial interval

include("Utilities/plots.jl")
include("Utilities/fileio.jl")


function main()

    # Select the target object type
    id = "Cube"

    # Read from corresponding directory
    dir = string("RealFlowData/", id, "/")
    fnames = readdir(dir)
    fnames = filter_unwanted_filenames(fnames)
    sort!(fnames, lt=trial_order)

    # Read all particle files
    dir = string("BulletData/", id, "/Intermediate/")
    all_files = filter_unwanted_filenames(readdir(dir)) 
    sort!(all_files, lt=trial_particle_order)

    #for i in eachindex(fnames)
    for i=0:9

        ela = 100*i
        
        for var=1:4

            trial_index = ela + var

            # Generate plot base
            plt1 = plot3d( 
                1,
                xlim = (-.5, .5),
                ylim = (-.5, .5),
                zlim = (0, 1),
                title = fnames[trial_index],
                legend = false,
                marker = 2,
                seriestype=:scatter
            )

            # Read ground truth trajectory
            head, tail = split(fnames[trial_index], '.')
            fname = join([head, "_observed.", tail])
            gt_file = string("RealFlowData/", id, "/", fname)
            print(string(gt_file, "\n"))
            ground_truth = CSV.read(gt_file, DataFrame)
            true_x = []
            true_y = []
            true_z = []

            # For each time step
            for x=1:30

                # Extend ground truth trajectory by one time step and add to plot
                true_x = [true_x; ground_truth[x, 1]]
                true_y = [true_y; ground_truth[x, 2]]
                true_z = [true_z; ground_truth[x, 3]]
                plot!(plt1, true_x, true_y, true_z, linewidth=3, linecolor=:red)

                # Extract the particle states at this time step
                particle_index = 30 * (trial_index-1) + x
                file = all_files[particle_index]
                print(string(file, "\n"))
                data = CSV.read(string(dir, file), DataFrame)
                tokens = split(file, "_")
                time_step = parse(Int64, replace(tokens[4], ".csv" => ""))
            
                # For each particle
                for i=1:20    

                    particle = data[data.particle .== i, :]

                    # elasticity = data[data.particle .== i, 2]
                    # weight = data[data.particle .== i, 3]
                    # @df plot!(plt1, particle[:, 5:7])

                    # Concatenate time steps to generate particle trajectory and add to plot
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
                if x % 5 == 0
                    display(plt1)
                    savefig(string("Plots/", id, "/", tokens[2], "_", tokens[3], "_", x))
                end
            end
        end
    end
end

main()
