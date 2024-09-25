# Data analysis
# Synthesizes and displays particle filter state at each time step
# Generates 3D plots showing particle trajectories vs ground truth
#
# QUESTION: Are the 'identities' of particles consistent across time, i.e. do they survive resampling?
# DONE: Allow more flexible selection of elasticity / trial interval
#
# TODO: Streamline display of dataframes (particle filter states) at each timestep
# TODO: Set the alpha / intensity to reflect the log weight of each particle

using Plots
using ZipFile

include("../Utilities/fileio.jl")
include("../Utilities/plots.jl")


# Model variation and experiment
model_id = "Modelv5"
target_id = "Sphere"
noise_id = "PosVar05"
expt_id = "BulletxBullet"
data_id = string(model_id, "/", target_id, "/", noise_id, "/", expt_id)

# Ground truth source file
gt_source = "Bullet"
realflow_expt = "Test"
gt_source == "Bullet" ? gt_dir = string("Tests/BulletStimulus/Data/", target_id, "/") : gt_dir = string("Data/RealFlowData/", realflow_expt, "/")
println(gt_dir)

# Plotting variables
interactive = false
plot_interval = 5


#########
# PLOTS #
#########

function make_plot_directories()
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
end

# Plot the trajectories 
function plot_trajectories(selected_files, expt_dir)

    # Plot parameters
    if interactive
        pyplot()
    else 
        make_plot_directories()
    end

    for (i, file) in enumerate(selected_files)

        println(file)

        # Read ground truth file to dataframe      
        gt_file = string(gt_dir, file)   
        ground_truth = CSV.read(gt_file, DataFrame)

        # Generate plot base
        tokens = split(file, "_")
        title = string("Stimulus: ", tokens[2], " ", tokens[3], "\n", data_id)
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

        # Procedurally generate plot
        true_x = []
        true_y = []
        true_z = []
        num_timesteps = size(ground_truth)[1]
        for t = 1:num_timesteps
            println(string("Timestep: ", t))

            # Extend ground truth trajectory by one time step and add to plot
            true_x = [true_x; ground_truth[t, 1]]
            true_y = [true_y; ground_truth[t, 2]]
            true_z = [true_z; ground_truth[t, 3]]
            plot!(plt, true_x, true_y, true_z, linewidth=3, linecolor=:red)

            # Every fifth time step, display the plot and (if static) save as .png
            if t % plot_interval == 0

                # Update title
                title!(string(title, "\nTimestep: ", t, " / ", num_timesteps))

                display(plt)

                tokens = split(file, "_")
                if !interactive
                    savefig(string("Analysis/Plots/Modelv5/Sphere/PosVar05/BulletxBullet/Trajectories/", tokens[2], "_", tokens[3], "_", t))
                end
            end
        end
    end
end


########
# MAIN #
########

function main()

    # Load and sort ground truth trajectory files
    gt_files = filter(contains("observed"), readdir(gt_dir))
    sort!(gt_files, lt=trial_order)

    plot_trajectories(gt_files, gt_dir)
    
end

# function run_with_profiler()
# # @profview main()
# @profile main()
# end
main()