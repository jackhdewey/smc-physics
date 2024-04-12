# Data analysis
# Generates 3D plots showing particle trajectories vs ground truth, with ground truth plotted in red
#
# DONE: Save plots in a folder
# DONE: Allow more flexible selection of elasticity / trial interval
#
# TODO: Set the alpha / intensity to reflect the log weight

include("Utilities/plots.jl")
include("Utilities/fileio.jl")


function main()

    # Select the target object type
    stimulus_id = "Cube"
    output_id = "SpherexCube"

    # Trial parameters
    num_timesteps = 30
    num_particles = 20

    # Other parameters
    elasticities = 10
    num_trials = 4
    plot_interval = 5

    # Read from corresponding directory
    dir = string("RealFlowData/", stimulus_id, "/")
    gt_files = filter_unwanted_filenames(readdir(dir))
    sort!(gt_files, lt=trial_order)

    # Read corresponding particle files
    dir = string("BulletData/", output_id, "/Intermediate/")
    particle_files = filter_unwanted_filenames(readdir(dir)) 
    sort!(particle_files, lt=trial_particle_order)

    #for i in eachindex(gt_files)

    # For each elasticity
    for ela=1:elasticities

        start_trial = 100 * (ela-1)
        
        # For first four trials
        for var=1:num_trials

            trial_index = start_trial + var

            # Generate plot base
            plt1 = plot3d( 
                1,
                xlim = (-.5, .5),
                ylim = (-.5, .5),
                zlim = (0, 1),
                title = gt_files[trial_index],
                legend = false,
                marker = 2,
                seriestype=:scatter
            )

            # Index into correct ground truth file
            head, tail = split(gt_files[trial_index], '.')
            fname = join([head, "_observed.", tail])
            gt_file = string("RealFlowData/", stimulus_id, "/", fname)
            print(string(gt_file, "\n"))

            # Read ground truth trajectory
            ground_truth = CSV.read(gt_file, DataFrame)
            true_x = []
            true_y = []
            true_z = []

            # For each time step
            for x=1:num_timesteps

                # Extend ground truth trajectory by one time step and add to plot
                true_x = [true_x; ground_truth[x, 1]]
                true_y = [true_y; ground_truth[x, 2]]
                true_z = [true_z; ground_truth[x, 3]]
                plot!(plt1, true_x, true_y, true_z, linewidth=3, linecolor=:red)

                # Index into correct particle filter file
                particle_index = 30 * (trial_index-1) + x
                file = particle_files[particle_index]

                # Extract the particle states at this time step
                data = CSV.read(string(dir, file), DataFrame)
                tokens = split(file, "_")
                time_step = parse(Int64, replace(tokens[4], ".csv" => ""))
            
                # For each particle
                for i=1:num_particles

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

                # Plot the particles at every fifth time step
                if x % plot_interval == 0
                    display(plt1)
                    savefig(string("Plots/", output_id, "/", tokens[2], "_", tokens[3], "_", x))
                end

            end
        end
    end
end

main()
