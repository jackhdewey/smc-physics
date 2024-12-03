# Data analysis
# Synthesizes and displays particle filter state at each time step
# Generates 3D plots showing particle trajectories vs ground truth
#
# DONE: Allow more flexible selection of elasticity / trial interval
# TODO: Streamline display of dataframes (particle filter states) at each timestep
# TODO: Set the alpha / intensity to reflect the log weight of each particle
#
# QUESTION: Are the 'identities' of particles consistent across time, i.e. do they survive resampling?

using Plots
using ZipFile

include("../args.jl")
include("../Utilities/fileio.jl")
include("../Utilities/plots.jl")


# Data source
gt_source = "Bullet"
obj_type = "Cube"

# Model parameters
model_id = "Modelv5"
target_id = "Cube"
noise_id = "PosVar075"

# Inference parameters
num_particles = 20

# Experiment
expt_id = "BulletTest"

# Plotting variables
interactive = false
plot_interval = 5

# Output filepath
data_path = string(model_id, "/", target_id, "/", noise_id, "/", expt_id)

# Displays (prints) the given particle filter states for the specified time steps
function display_data_frames(r, particle_indices, interval)

    # For the files starting at the selected indices
    for particle_index in particle_indices

        # For an interval of particle filter time steps
        for t=interval[1]:interval[2]
            file = r.files[particle_index + t - 1]
            data = CSV.File(read(file)) |> DataFrame

            # Select only one frame for each particle
            time_step = data[data.frame.==1, :]

            # For each particle, check if its elasticity is unique
            unique = Dict()
            for i=1:num_particles
                if !haskey(unique, time_step[i, 2])
                    unique[time_step[i, 2]] = 1
                else   
                    unique[time_step[i, 2]] += 1
                end
            end
            println(unique)
        end
    end

end

# Plot the trajectories 
function plot_trajectories(gt_files, particle_indices, r, expt_dir)

    # Plot parameters
    if interactive
        pyplot()
    end

    for i in eachindex(gt_files)

        # Read ground truth file to dataframe
        gt_file = string(gt_dir, gt_files[i])
        ground_truth = CSV.read(gt_file, DataFrame)

        # Generate plot base
        tokens = split(gt_files[i], "_")
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

            # Extract and sort particle files
            files = map((file) -> file.name, r.files)
            files = filter(contains(".csv"), files)
            sort!(files, lt=trial_particle_order)

            # Index into correct file and read into dataframe
            d_file = files[particle_indices[i] + t - 1]
            t_file = filter((file) -> file.name == d_file, r.files)[1]
            data = CSV.File(read(t_file)) |> DataFrame

            # For each particle
            for p = 1:num_particles
                particle = data[data.particle.==p, :]
                println(particle)

                # elasticity = data[data.particle .== i, 2]
                # weight = data[data.particle .== i, 3]
                # @df plot!(plt1, particle[:, 5:7])

                # Read particle trajectory up to current time step
                x_trajectory = []
                y_trajectory = []
                z_trajectory = []
                for row = 1:t
                    x_trajectory = [x_trajectory; particle[row, 5]]
                    y_trajectory = [y_trajectory; particle[row, 6]]
                    z_trajectory = [z_trajectory; particle[row, 7]]
                end

                # Add particle trajectory to plot
                plot!(plt, x_trajectory, y_trajectory, z_trajectory)

            end

            # Every fifth time step, display the plot and (if static) save as .png
            if t % plot_interval == 0

                # Update title
                title!(string(title, "\nTimestep: ", t, " / ", num_timesteps))

                display(plt)

                tokens = split(t_file.name, "_")
                if !interactive
                    savefig(string(expt_dir, "/Trajectories/", tokens[1], "_", tokens[2], "_", t))
                end
            end
        end
        #particle_index += num_timesteps
    end
end

# Check whether the output filepath already exists, if not make directories
function make_output_directories()

    plot_dir = string("Analysis/Plots/", model_id, "/")
    if !isdir(plot_dir)
        mkdir(plot_dir)
    end
    target_dir = string(plot_dir, target_id, "/")
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
    data_dir = string(expt_dir, "Trajectories/")
    if !isdir(data_dir)
        mkdir(data_dir)
    end

end


function main()

    # Load and sort ground truth trajectory files
    gt_source == "Bullet" ? 
        gt_dir = string("Tests/BulletStimulus/Data/", obj_type, "/") : 
        gt_dir = string("Data/RealFlowData/", expt_id, "/")
    gt_files = filter(contains("observed"), readdir(gt_dir))
    sort!(gt_files, lt=trial_order)

    # Load and sort intermediate particle filter state files
    dir = string("Data/BulletData/", data_path, "/Intermediate/")
    r = ZipFile.Reader(string(dir, "particles.zip"))
    files = map((file) -> file.name, r.files)
    particle_files = filter(contains(".csv"), files)
    sort!(particle_files, lt=trial_particle_order)

    make_output_directories()

    # Filter to trials we want to plot
    #=
    #gt_files = filter((file) -> occursin("Ela9_Var1_", file), gt_files)

    #_, error_trials = process_individual_stimuli_sim(1, target_id)

    error_trials = [
        "Sphere_Ela0_Var4", "Sphere_Ela0_Var6", "Sphere_Ela0_Var14", 

        "Sphere_Ela1_Var10", "Sphere_Ela1_Var12", "Sphere_Ela1_Var2", "Sphere_Ela1_Var4", "Sphere_Ela1_Var6", 

        "Sphere_Ela2_Var12", "Sphere_Ela2_Var13", "Sphere_Ela2_Var6", 
        "Sphere_Ela2_Var7", 

        "Sphere_Ela3_Var1", "Sphere_Ela3_Var13", "Sphere_Ela3_Var14", "Sphere_Ela3_Var5", "Sphere_Ela3_Var8", 

        "Sphere_Ela4_Var1", "Sphere_Ela4_Var14", "Sphere_Ela4_Var15", "Sphere_Ela4_Var2", "Sphere_Ela4_Var5", 

        "Sphere_Ela5_Var10", "Sphere_Ela5_Var11", "Sphere_Ela5_Var2", "Sphere_Ela5_Var5", 

        "Sphere_Ela6_Var10", "Sphere_Ela6_Var11", "Sphere_Ela6_Var12", "Sphere_Ela6_Var13", "Sphere_Ela6_Var14", "Sphere_Ela6_Var2", "Sphere_Ela6_Var4", "Sphere_Ela6_Var9", 

        "Sphere_Ela7_Var10", "Sphere_Ela7_Var2", "Sphere_Ela7_Var6"
    ]
    sort!(error_trials, lt=trial_order)

    # Locate the indices of the corresponding RealFlow files
    trial_indices = []
    for file in error_trials
        tokens = split(file, "_")
        key = string(tokens[2], "_", tokens[3], "_")
        trial_index = findfirst(file -> occursin(key, file), gt_files)
        push!(trial_indices, trial_index)
    end

    # Find the indices of the corresponding particle filter files
    particle_indices = []
    particle_index = 1
    for i = 1:maximum(trial_indices)
        if i in trial_indices
            push!(particle_indices, particle_index)
        end
        gt_file = string(gt_dir, gt_files[i])
        particle_index += size(CSV.read(gt_file, DataFrame))[1]
    end
 
    # Check to confirm corresponding file names
    for i in eachindex(trial_indices)
        println(error_trials[i])
        println(files[particle_indices[i]])
    end
    =#

    # Locate the indices of the intermediate particle filter state files for each gt_file
    particle_index = 1
    particle_indices = []
    for file in gt_files
        
        push!(particle_indices, particle_index)

        gt_file = string(gt_dir, file)
        frames = size(CSV.read(gt_file, DataFrame))[1]
        println(frames)
        particle_index += frames

    end

    interval=[1, 10]
    #display_data_frames(r, particle_indices, interval)
    plot_trajectories(gt_files, particle_indices, r, expt_dir)
    
end

# function run_with_profiler()
# # @profview main()
# @profile main()
# end
main()
