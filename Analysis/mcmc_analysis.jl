# Data analysis
# Synthesizes and displays MCMC state for each time run
# Generates 3D plots showing simulated trajectories vs ground truth

using Plots
using ZipFile

include("../args.jl")
include("../Utilities/fileio.jl")
include("../Utilities/plots.jl")


# Ground truth source file
gt_source = "Bullet"

# Model variation 
model_id = "Modelv5"
target_id = "Sphere"
noise_id = "PosVar05"

expt_id = "BulletTest"

# Output filepath
data_path = string(model_id, "/", target_id, "/", noise_id, "/", expt_id)

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
    alg_dir = string(noise_dir, "MCMC/")
    if !isdir(alg_dir)
        mkdir(alg_dir)
    end
    expt_dir = string(alg_dir, expt_id, "/")
    if !isdir(expt_dir)
        mkdir(expt_dir)
    end
    output_dir = string(expt_dir, "Trajectories/")
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    return output_dir
end

function plot_ground_truth(gt_dir, file, output_dir)

    # Read ground truth file to dataframe      
    gt_file = string(gt_dir, file)   
    ground_truth = CSV.read(gt_file, DataFrame)

    # Generate plot base
    tokens = split(file, "_")
    title = string("Stimulus: ", tokens[2], " ", tokens[3], "\n", output_dir)
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
                savefig(string(output_dir, tokens[2], "_", tokens[3], "_", t))
            end

        end
    end

    return plt

end

# Plot the trajectories 
function plot_trajectories(plt, r)

    # Extract and sort particle files
    files = map((file) -> file.name, r.files)
    mcmc_files = filter(contains(".csv"), files)
    sort!(mcmc_files, lt=trial_order)

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
                savefig(string(output_dir, tokens[2], "_", tokens[3], "_", t))
            end
        end
    end

end


########
# MAIN #
########

function main()

    # Load and sort ground truth trajectory files
    gt_source == "Bullet" ? 
        gt_dir = string("Tests/BulletStimulus/Data/", target_id, "/") : 
        gt_dir = string("Data/RealFlowData/", realflow_expt, "/")
    println(gt_dir)
    gt_files = filter(contains("observed"), readdir(gt_dir))
    sort!(gt_files, lt=trial_order)

    inference_dir = string("Data/BulletData/", data_path, "Trajectories")
    r = ZipFile.Reader(string(dir, "inferences.zip"))

    # Plot parameters
    if interactive
        pyplot()
    else 
        output_dir = make_plot_directories()
    end

    for file in gt_files
        plt = plot_ground_truth(gt_dir, file, output_dir)
        #plot_trajectories(plt, r)
    end
    
end

# function run_with_profiler()
# # @profview main()
# @profile main()
# end
main()