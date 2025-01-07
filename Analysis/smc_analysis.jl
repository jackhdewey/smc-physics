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


# Plotting variables
interactive = false
plot_interval = 5

# Displays (prints) the given particle filter states for the specified time steps
function display_data_frames(reader, particle_indices, interval)

    # For the files starting at the selected indices
    for particle_index in particle_indices

        # For an interval of particle filter time steps
        for t=interval[1]:interval[2]
            file = reader.files[particle_index + t - 1]
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
function plot_trajectories(gt_files, particle_indices, reader, plot_path)

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
            files = map((file) -> file.name, reader.files)
            files = filter(contains(".csv"), files)
            sort!(files, lt=trial_particle_order)

            # Index into correct particle file and read into dataframe
            d_file = files[particle_indices[i] + t - 1]
            t_file = filter((file) -> file.name == d_file, reader.files)[1]
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
                    savefig(string(plot_path, tokens[1], "_", tokens[2], "_", t))
                end
            end
        end
        
        # Update index to correct particle filter file
        particle_index += num_timesteps

    end
end


function main()

    args = Args()

    # Load and sort ground truth trajectory files
    args.gt_source == "Bullet" ? 
        gt_dir = string("Tests/BulletStimulus/Data/", args.gt_shape, "/") : 
        gt_dir = string("Data/RealFlowData/", args.expt_id, "/")
    gt_files = filter(contains("observed"), readdir(gt_dir))
    sort!(gt_files, lt=trial_order)

    noise_id = generate_noise_id(args)
    inference_id = generate_inference_param_id(args)

    # Load and sort intermediate particle filter state files
    data_path = joinpath("Analysis", args.expt_id, args.model_id, args.target_id, noise_id, args.algorithm, "Data")
    files = readdir(data_path)
    
    #reader = ZipFile.Reader(string(dir, "particles.zip"))
    #files = map((file) -> file.name, reader.files)

    particle_files = filter(contains(".csv"), files)
    sort!(particle_files, lt=trial_particle_order)

    # For each gt_file, store the index of the first corresponding intermediate particle filter state file
    particle_index = 1
    particle_indices = []
    for file in gt_files
        
        push!(particle_indices, particle_index)
    
        gt_file = string(gt_dir, file)
        num_frames = size(CSV.read(gt_file, DataFrame))[1]
        println(num_frames)
        particle_index += frames
    
    end

    plot_path = generate_plot_path(args.expt_id, args.model_id, args.target_id, noise_id, inference_id, "Trajectories")
    plot_trajectories(gt_files, particle_indices, reader, plot_path)

    #interval=[1, 10]
    #display_data_frames(reader, particle_indices, interval)

    #make_output_directories()

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
    
end

# function run_with_profiler()
# # @profview main()
# @profile main()
# end
main()
