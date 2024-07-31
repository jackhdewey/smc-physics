# Data analysis
# Generates 3D plots showing particle trajectories vs ground truth, with ground truth plotted in red
#
# DONE: Save plots in a folder
# DONE: Allow more flexible selection of elasticity / trial interval
#
# TODO: Set the alpha / intensity to reflect the log weight of each particle

include("../Utilities/fileio.jl")
include("../Utilities/plots.jl")

using Plots
# using PyPlot
pyplot()

function main()

    # Select the model variation and experiment
    model_id = "Modelv5"
    target_id = "Sphere"
    noise_id = "PosVar075"
    expt_id = "Test"
    output_id = string(model_id, "/", target_id, "/", noise_id, "/", expt_id)

    # Plot parameters
    interactive = true
    if interactive
        pyplot()
    end
    plot_interval = 5

    # Pull ground truth trajectory files from directory
    dir = string("Data/RealFlowData/", expt_id, "/")
    gt_files = filter(contains("observed"), readdir(dir))
    sort!(gt_files, lt=trial_order)

    # Pull intermediate particle filter state files from directory
    dir = string("Data/BulletData/", output_id, "/Intermediate/")
    particle_files = filter_unwanted_filenames(readdir(dir))
    sort!(particle_files, lt=trial_particle_order)

    data = CSV.read(string(dir, particle_files[1]), DataFrame)
    num_particles = size(data)[1]

    # For each trial
    total_particle_index = 0
    for i in eachindex(gt_files)

        plotp = false


        if occursin("Ela9_Var1_", gt_files[i])
            plotp = true
        end
        # Read ground truth file to dataframe
        gt_file = string("Data/RealFlowData/", expt_id, "/", gt_files[i])
        ground_truth = CSV.read(gt_file, DataFrame)

        # Generate plot base
        tokens = split(gt_files[i], "_")
        if plotp
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
        )
        end

        # Procedurally generate plots for each timestep
        num_timesteps = size(ground_truth)[1]
        # true_x = zeros(num_timesteps)
        # true_y = zeros(num_timesteps)
        # true_z = zeros(num_timesteps)
        true_x = []
        true_y = []
        true_z = []
        for t = 1:num_timesteps

            # Extend ground truth trajectory by one time step and add to plot
            true_x = [true_x; ground_truth[t, 1]]
            true_y = [true_y; ground_truth[t, 2]]
            true_z = [true_z; ground_truth[t, 3]]
            # true_x[t] = ground_truth[t, 1]
            # true_y[t] = ground_truth[t, 2]
            # true_z[t] = ground_truth[t, 3]
            if plotp
            plot!(plt, true_x, true_y, true_z, linewidth=3, linecolor=:red)
            end

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

                # Generate particle trajectory up to current time step and add to plot
                x_trajectory = []
                y_trajectory = []
                z_trajectory = []
                # x_trajectory[row] = particle[row, 5]
                # y_trajectory[row] = particle[row, 6]
                # z_trajectory[row] = particle[row, 7]
                for row = 1:t
                    x_trajectory = [x_trajectory; particle[row, 5]]
                    y_trajectory = [y_trajectory; particle[row, 6]]
                    z_trajectory = [z_trajectory; particle[row, 7]]
                end
                if plotp
                plot!(plt, x_trajectory, y_trajectory, z_trajectory)
                #PyPlot.text3D(0, 0, 0, "test")
                #annotate!(x_trajectory[end], y_trajectory[end], z_trajectory[end], text(t, 10))
                # PyPlot.text3D(x_trajectory[end], y_trajectory[end], z_trajectory[end], "1")
                end

            end

            # Display the plot every fifth time step and save if static
            tokens = split(file, "_")
            if plotp
            if t % plot_interval == 0
                display(plt)
                if !interactive
                    savefig(string("Analysis/Plots/", output_id, "/", tokens[2], "_", tokens[3], "_", t))
                end
            end
        end
        end
        total_particle_index += num_timesteps
    end
end

# function run_with_profiler()
# # @profview main()
# @profile main()
# end
main()
