# Data analysis
# Generates 3D plots showing particle trajectories vs ground truth, with ground truth plotted in red
#
# DONE: Allow more flexible selection of elasticity / trial interval
#
# TODO: 
# TODO: Set the alpha / intensity to reflect the log weight of each particle

include("../Utilities/fileio.jl")
include("../Utilities/plots.jl")

using Plots
using ZipFile

function main()

    # Select the model variation and experiment
    model_id = "Modelv5"
    target_id = "Cube"
    noise_id = "PosVar075"
    expt_id = "Test"
    data_id = string(model_id, "/", target_id, "/", noise_id, "/", expt_id)

    num_particles = 20

    # Plot parameters
    interactive = false
    if interactive
        pyplot()
    end
    plot_interval = 5

    # Pull ground truth trajectory files from directory
    dir = string("Data/RealFlowData/", expt_id, "/")
    gt_files = filter(contains("observed"), readdir(dir))
    sort!(gt_files, lt=trial_order)

    # Pull intermediate particle filter state files from directory
    dir = string("Data/BulletData/", data_id, "/Intermediate/")
    r = ZipFile.Reader(string(dir, "particles.zip"))

    #particle_files = map(file -> file.name, r.files)
    #sort!(particle_files, lt=trial_particle_order)

    model_dir = string("Analysis/Plots/", model_id, "/")
    if !isdir(model_dir)
        mkdir(model_dir)
    end
    target_dir = string(model_dir, target_id, "/")
    if !isdir(target_dir)
        mkdir(target_dir)
    end
    noise_dir = string(target_dir, noise_id, "/")
    if !isdir(noise_dir)
        mkdir(noise_dir)
    end
    expt_dir = string(noise_dir, expt_id, "/")
    if !isdir(expt_dir)
        mkdir(expt_dir)
    end

    # Filter which trials we want to plot
    trial_index = findall(file -> occursin("Ela9_Var1_", file), gt_files)
    println(trial_index[1])
    println(gt_files[trial_index])

    particle_index = 1
    for i = 1:trial_index[1]-1
        gt_file = string("Data/RealFlowData/", expt_id, "/", gt_files[i])
        particle_index += size(CSV.read(gt_file, DataFrame))[1]
    end
    println(particle_index)

    # For each trial
    gt_files = filter((file) -> occursin("Ela9_Var1_", file), gt_files)
    for i in eachindex(gt_files)

        tokens = split(gt_files[i], "_")
        title = string("Stimulus: ", tokens[2], " ", tokens[3], "\n", data_id)

        # Generate plot base
        plt = plot3d(
                1,
                xlim=(-0.5, 0.5),
                ylim=(-0.5, 0.5),
                zlim=(0, 1),
                title=title,
                legend=false,
                marker=2,
                seriestype=:scatter,
                size=(1200, 800),
                gridlinewidth=8
        )

        # Read ground truth file to dataframe
        gt_file = string("Data/RealFlowData/", expt_id, "/", gt_files[i])
        ground_truth = CSV.read(gt_file, DataFrame)

        # Procedurally generate plots for each timestep
        num_timesteps = size(ground_truth)[1]
        true_x = []
        true_y = []
        true_z = []
        for t = 1:num_timesteps
            println(string("Timestep: ", t))

            # Extend ground truth trajectory by one time step and add to plot
            true_x = [true_x; ground_truth[t, 1]]
            true_y = [true_y; ground_truth[t, 2]]
            true_z = [true_z; ground_truth[t, 3]]
            
            plot!(plt, true_x, true_y, true_z, linewidth=3, linecolor=:red)

            # Index into correct particle filter file
            file = r.files[particle_index + t - 1]
            println(file.name)
            data = CSV.File(read(file)) |> DataFrame

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
                for row = 1:t
                    x_trajectory = [x_trajectory; particle[row, 5]]
                    y_trajectory = [y_trajectory; particle[row, 6]]
                    z_trajectory = [z_trajectory; particle[row, 7]]
                end
                title!(string(title, "\nTimestep: ", t, " / ", num_timesteps))
                plot!(plt, x_trajectory, y_trajectory, z_trajectory)

            end

            # Display the plot every fifth time step and save if static
            tokens = split(file.name, "_")
            if t % plot_interval == 0
                display(plt)
                if !interactive
                    savefig(string(expt_dir, "/", tokens[2], "_", tokens[3], "_", t))
                end
            end
        end
        particle_index += num_timesteps
    end
    
end

# function run_with_profiler()
# # @profview main()
# @profile main()
# end
main()
