# Runs analysis scripts and generates plots 

include("../args.jl")
include("plots.jl")
include("utils.jl")


# Generate paths to relevant data
function get_paths(args)

    # Load and sort ground truth trajectory files
    if contains(args.expt_id, "Bullet")  
        bullet_shape = split(args.expt_id, "_")[2] 
        gt_path = string("Tests/BulletStimulus/Data/", bullet_shape, "/") 
    else
        gt_path = string("Data/RealFlowData/", args.expt_id, "/") 
    end

    noise_id = generate_noise_id(args)
    inference_id = generate_inference_param_id(args)
    trial_path = joinpath("Analysis", args.expt_id, args.model_id, args.target_id, noise_id, args.algorithm, inference_id)
    plots_path = generate_plot_path(args, noise_id, inference_id)

    return gt_path, trial_path, plots_path
end

# Read ground truth trajectory data and intermediate particle filter data
function get_files(gt_path, trial_path)

    # Load and sort ground truth trajectory files
    gt_files = filter(contains("observed"), readdir(gt_path))
    sort!(gt_files, lt=trial_order)

    # Generate path(s) to data
    data_path = joinpath(trial_path, "Data")
    particle_data_path = joinpath(data_path, "intermediate")

    # Load and sort intermediate particle filter state files
    particle_files = readdir(particle_data_path)
    particle_files = filter(contains(".csv"), particle_files)
    #reader = ZipFile.Reader(string(data_path, "particles.zip"))
    #particle_files = map((file) -> file.name, reader.files)
    sort!(particle_files, lt=trial_particle_order)

    # Load inference output files
    output_files = readdir(joinpath(data_path, "inferences"))
    
    return gt_files, particle_files, output_files

end

# Identify trials with high model error (vs. human)
# Evaluate and plot model performance at inferring elasticity relative to ground truth
# Evaluate and plot model-inferred vs human elasticity judgments
function analyze_inference(sim_data, plots_path, args)

    # Extract high-error trials
    sim_data, error_trials = process_individual_stimuli_sim(sim_data)
    println("High Error Trials: ", error_trials)

    # Generate output filepath for plots
    plots_path = joinpath(plots_path, "Judgments")
    println("Plots Path: ", plots_path)

    # Set plot marker shape
    if args.target_id == "Cube"
        marker_shape = :square
    else
        marker_shape = :circle
    end

    # Plot model-inferred elasticity against ground truth
    plot_inferences_vs_gt("sim", sim_data, args.expt_id, args.target_id, marker_shape, plots_path)

    # RealFlow trials only
    if occursin("Exp", args.expt_id)

        # Read human data into data frame
        human_data = read_subject_data(args.expt_id)
        human_data = process_individual_stimuli_human(human_data)

        # Plot human judgments against ground truth
        plot_inferences_vs_gt("human", human_data, args.expt_id, args.target_id, marker_shape, plots_path)

        # Plot model judgments against human
        plot_inferences_sim_vs_human(sim_data, human_data, args.expt_id, args.target_id, marker_shape, plots_path)

    end

end

# Generates 3D plots showing particle trajectories and ground truth, i.e. displaying particle filter state at each time step
function analyze_trajectories(particle_files, gt_files, plots_path, args)

    plot_path = joinpath(plots_path, "Trajectories")

    # Store the index of the first intermediate particle filter state file corresponding to each gt_file
    particle_indices = []
    particle_index = 1
    for file in gt_files
        
        push!(particle_indices, particle_index)
        gt_file = string(gt_dir, file)
        num_frames = size(CSV.read(gt_file, DataFrame))[1]
        particle_index += num_frames
    
    end

    plot_trajectories(gt_dir, gt_files, data_path, particle_files, particle_indices, args.num_particles, plots_path)

    #interval=[1, 10]
    #display_data_frames(reader, particle_indices, interval)

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

function main()

    args = Args()
    gt_path, trial_path, plots_path = get_paths(args)
    gt_files, particle_files, output_files = get_files(gt_path, trial_path)
    output_data = read_simulation_data(output_files)
    
    #analyze_trajectories(particle_files, gt_files, plots_path, args)
    generate_violin_plot(args.expt_id, args.model_id, args.param_id, args.num_particles, plots_path)
    #analyze_inference(output_data, plots_path, args)

end

main()
