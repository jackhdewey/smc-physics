# Generates 3D plots showing particle trajectories and ground truth, i.e. displaying particle filter state at each time step

# QUESTION: Are the 'identities' of particles consistent across time, i.e. do they survive resampling?

# TODO: Set the color based on the elasticity of each particle
# TODO: Set the alpha / intensity to reflect the log weight of each particle
# TODO: Streamline display of dataframes (particle filter states) at each timestep

using Plots
using ZipFile

include("../args.jl")
include("../Utilities/fileio.jl")
include("../Utilities/plots.jl")
include("utils.jl")


# Plotting variables
interactive = false
plot_interval = 5

function main()

    args = Args()

    # Load and sort ground truth trajectory files
    if contains(args.expt_id, "Bullet")  
        bullet_shape = split(args.expt_id, "_")[2] 
        gt_dir = string("Tests/BulletStimulus/Data/", bullet_shape, "/") 
    else
        gt_dir = string("Data/RealFlowData/", args.expt_id, "/") 
    end
    gt_files = filter(contains("observed"), readdir(gt_dir))
    sort!(gt_files, lt=trial_order)

    # Load and sort intermediate particle filter state files
    noise_id = generate_noise_id(args)
    inference_id = generate_inference_param_id(args)
    data_path = joinpath("Analysis", args.expt_id, args.model_id, args.target_id, noise_id, args.algorithm, inference_id, "Data", "intermediate")
    #reader = ZipFile.Reader(string(dir, "particles.zip"))
    #files = map((file) -> file.name, reader.files)
    particle_files = readdir(data_path)
    particle_files = filter(contains(".csv"), particle_files)
    sort!(particle_files, lt=trial_particle_order)

    # For each gt_file, store the index of the first corresponding intermediate particle filter state file
    particle_index = 1
    particle_indices = []
    for file in gt_files
        
        push!(particle_indices, particle_index)
    
        gt_file = string(gt_dir, file)
        num_frames = size(CSV.read(gt_file, DataFrame))[1]
        particle_index += num_frames
    
    end

    plot_path = generate_plot_path(args.expt_id, args.model_id, args.target_id, noise_id, inference_id, "Trajectories")
    plot_trajectories(gt_dir, gt_files, data_path, particle_files, particle_indices, args.num_particles, plot_path)

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

# function run_with_profiler()
# # @profview main()
# @profile main()
# end
main()
