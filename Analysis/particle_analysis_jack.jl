# Data analysis
# Generates 3D plots showing particle trajectories vs ground truth
#
# DONE: Allow more flexible selection of elasticity / trial interval
# TODO: Streamline display of dataframes (particle filter states) at each timestep to assess behavior
# TODO: Set the alpha / intensity to reflect the log weight of each particle

using Plots
using ZipFile

include("../Utilities/fileio.jl")
include("../Utilities/plots.jl")

# Select the model variation and experiment
model_id = "Modelv5"
target_id = "Sphere"
noise_id = "PosVar075"
expt_id = "Exp1"
data_id = string(model_id, "/", target_id, "/", noise_id, "/", expt_id)

num_particles = 20

# Display data frames for the specified interval of particle filter states
function display_data_frames(r, particle_indices, interval)

    for particle_index in particle_indices
        # For an interval of particle filter time steps
        for t=interval[1]:interval[2]
            file = r.files[particle_index + t - 1]
            data = CSV.File(read(file)) |> DataFrame

            # For the first time step
            time_step = data[data.frame.==1, :]
            unique = Dict()

            # For each particle
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
    interactive = false
    if interactive
        pyplot()
    end
    plot_interval = 5

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

        # Procedurally generate plots by adding trajectories at each time step
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
            files = map((file) -> file.name, r.files)
            files = filter(contains(".csv"), files)
            sort!(files, lt=trial_particle_order)
            d_file = files[particle_indices[i] + t - 1]
            println(d_file)
            t_file = filter((file) -> file.name == d_file, r.files)[1]
            println(t_file.name)
            data = CSV.File(read(t_file)) |> DataFrame

            # For each particle
            for p = 1:num_particles
                particle = data[data.particle.==p, :]
                println(particle)

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

            # Every fifth time step, display the plot and (if static) save as .png
            tokens = split(t_file.name, "_")
            if t % plot_interval == 0
                display(plt)
                if !interactive
                    savefig(string(expt_dir, "/Trajectories/", tokens[2], "_", tokens[3], "_", t))
                end
            end
        end
        #particle_index += num_timesteps
    end
end

function main()

    # Pull ground truth trajectory files from directory
    dir = string("Data/RealFlowData/", expt_id, "/")
    gt_files = filter(contains("observed"), readdir(dir))
    sort!(gt_files, lt=trial_order)

    # Pull intermediate particle filter state files from directory
    dir = string("Data/BulletData/", data_id, "/Intermediate/")
    r = ZipFile.Reader(string(dir, "particles.zip"))
    files = map((file) -> file.name, r.files)
    files = filter(contains(".csv"), files)
    sort!(files, lt=trial_particle_order)

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

    # Filter down to trials we want to plot

    # For each trial
    #gt_files = filter((file) -> occursin("Ela9_Var1_", file), gt_files)

    #_, error_trials = process_individual_stimuli_sim(1, target_id)
    error_trials = ["Sphere_Ela0_Var14", "Sphere_Ela0_Var4", "Sphere_Ela0_Var6", "Sphere_Ela1_Var10", "Sphere_Ela1_Var12", 
    "Sphere_Ela1_Var2", "Sphere_Ela1_Var4", "Sphere_Ela1_Var6", "Sphere_Ela2_Var12", "Sphere_Ela2_Var13", "Sphere_Ela2_Var6", 
    "Sphere_Ela2_Var7", "Sphere_Ela3_Var1", "Sphere_Ela3_Var13", "Sphere_Ela3_Var14", "Sphere_Ela3_Var5", "Sphere_Ela3_Var8", 
    "Sphere_Ela4_Var1", "Sphere_Ela4_Var14", "Sphere_Ela4_Var15", "Sphere_Ela4_Var2", "Sphere_Ela4_Var5", "Sphere_Ela5_Var10", "Sphere_Ela5_Var11", "Sphere_Ela5_Var2", "Sphere_Ela5_Var5", "Sphere_Ela6_Var10", "Sphere_Ela6_Var11", "Sphere_Ela6_Var12", "Sphere_Ela6_Var13", "Sphere_Ela6_Var14", "Sphere_Ela6_Var2", "Sphere_Ela6_Var4", "Sphere_Ela6_Var9", "Sphere_Ela7_Var10", "Sphere_Ela7_Var2", "Sphere_Ela7_Var6"]

    sort!(error_trials, lt=trial_order)
    trial_indices = []
    for file in error_trials
        tokens = split(file, "_")
        key = string(tokens[2], "_", tokens[3], "_")
        trial_index = findfirst(file -> occursin(key, file), gt_files)
        push!(trial_indices, trial_index)
    end

    #println("Error Trials: ", trial_indices)
    #=
    for i in trial_indices
        println(gt_files[i])
    end
    =#

    particle_indices = []
    particle_index = 1
    max_index = maximum(trial_indices)
    for i = 1:max_index
        if i in trial_indices
            push!(particle_indices, particle_index)
        end
        gt_file = string("Data/RealFlowData/", expt_id, "/", gt_files[i])
        particle_index += size(CSV.read(gt_file, DataFrame))[1]
    end
 
    for i in eachindex(trial_indices)
        println(error_trials[i])
        println(files[particle_indices[i]])
    end
    
    plot_trajectories(gt_files, particle_indices, r, expt_dir)

    #gt_file = string("Data/RealFlowData/", expt_id, "/", gt_files[1])
    #ground_truth = CSV.read(gt_file, DataFrame)

    #interval=[1, 10]
    #display_data_frames(r, particle_indices, interval)
    
end

# function run_with_profiler()
# # @profview main()
# @profile main()
# end
main()
