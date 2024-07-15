# Data analysis
# Generates 3D plots showing particle trajectories vs ground truth, with ground truth plotted in red
#
# DONE: Save plots in a folder
# DONE: Allow more flexible selection of elasticity / trial interval
#
# TODO: Compute error
# TODO: Set the alpha / intensity to reflect the log weight of each particle

include("../Utilities/fileio.jl")
include("../Utilities/plots.jl")

using Plots

pyplot()


function main()

    # Select the target object type
    model_id = "Modelv5"
    target_id = "Sphere"
    noise_id = "PosVar075"
    expt_id = "Test"
    output_id = string(model_id, "/", target_id, "/", noise_id, "/", expt_id)

    # Trial parameters
    num_particles = 20

    # Other parameters
    num_trials = 4
    plot_interval = 5

    # Pull ground truth (RealFlow) trajectory files from directory
    dir = string("Data/RealFlowData/", expt_id, "/")
    gt_files = filter(contains("observed"), readdir(dir))
    sort!(gt_files, lt=trial_order)
    # print(gt_files)

    # Pull intermediate particle filter state files from directory
    dir = string("Data/BulletData/", output_id, "/Intermediate/")
    particle_files = filter_unwanted_filenames(readdir(dir)) 
    sort!(particle_files, lt=trial_particle_order)

    # For each trajectory
    total_particle_index = 0

    for i in eachindex(gt_files[1:2])

    for i in eachindex(gt_files)[1:1]
        tokens = split(gt_files[i], "_")

        # Generate plot base
        plt = plot3d(
            1,
            xlim=(-0.5, 0.5),
            ylim=(-0.5, 0.5),
            zlim=(0, 1),
            title=string("Stimulus: ", tokens[2], " ", tokens[3], "\n", model_id),
            legend=false,
            marker=2,
            seriestype=:scatter,
            size=(1200, 800),
            gridlinewidth=8
            # aspect_ratio = 1
            # layout=layout
        )
        # p = PlotlyJSBackend().o

        # # Adjust the aspect ratio using PlotlyJS layout
        # layout = Layout(
        #     scene=attr(
        #         aspectmode="cube"  # Ensures equal aspect ratio
        #     )
        # )

        # relayout!(p, layout)
        gt_file = string("Data/RealFlowData/", expt_id, "/", gt_files[i])
        ground_truth = CSV.read(gt_file, DataFrame)

        # Generate the plot procedurally
        num_timesteps = size(ground_truth)[1]
        true_x = []
        true_y = []
        true_z = []
        for t = 1:num_timesteps
            # for t = 10:11

            # Extend ground truth trajectory by one time step and add to plot
            true_x = [true_x; ground_truth[t, 1]]
            true_y = [true_y; ground_truth[t, 2]]
            true_z = [true_z; ground_truth[t, 3]]
            plot!(plt, true_x, true_y, true_z, linewidth=3, linecolor=:red)

            # Index into correct particle filter file
            particle_index = total_particle_index + t
            file = particle_files[particle_index]
            data = CSV.read(string(dir, file), DataFrame)

            # For each particle
            for p = 1:num_particles

                particle = data[data.particle.==p, :]

                # elasticity = data[data.particle .== i, 2]
                # weight = data[data.particle .== i, 3]
                # @df plot!(plt1, particle[:, 5:7])

                # Generate complete particle trajectory up to current time step and add to plot
                x_trajectory = []
                y_trajectory = []
                z_trajectory = []
                for row = 1:t
                    x_trajectory = [x_trajectory; particle[row, 5]]
                    y_trajectory = [y_trajectory; particle[row, 6]]
                    z_trajectory = [z_trajectory; particle[row, 7]]
                end
                plot!(plt, x_trajectory, y_trajectory, z_trajectory)

            end

            # Save and display the plot every fifth time step
            tokens = split(file, "_")
            if t % plot_interval == 0
                # savefig(string("Analysis/Plots/", output_id, "/", tokens[2], "_", tokens[3], "_", t))
                display(plt)
            end

        end
        total_particle_index += num_timesteps
    end
end

main()
