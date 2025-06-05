# Utility functions for analysis scripts - mostly file I/O

using DataFrames
using DataFramesMeta
using CSV
using ZipFile

project_path = dirname(@__DIR__)


###############
#  READ DATA  #
###############

# Read all simulation data files at the specified directory 
function read_simulation_data(data_files, zip=true)
    
    all_data = []

    for file in data_files
        
        if zip 
            filename = file.name
        else
            filename = joinpath(simulation_folder, "inferences", file)
        end

        if filename == "inferences/"
            continue
        end

        # Extract ground truth elasticity from filename
        tokens = split(filename, ['/', '_', '.'])
        elasticity_string, variation_string = tokens[length(tokens)-2], tokens[length(tokens)-1]
        elasticity = parse(Int, elasticity_string[end]) * 0.1   # fix order of magnitude
        variation = parse(Int64, variation_string[4:end])

        # Try to read the file into a dataframe
        println("Reading file ", filename)
        try
            if zip
                data = CSV.File(read(file)) |> DataFrame
            else 
                data = CSV.File(read(filename)) |> DataFrame
            end
            data = insertcols(
                data,
                "filename" => file,
                "stimulusID" => target_id * "_" * elasticity_string * "_" * variation_string,
                "gtElasticity" => elasticity,
                "variation" => variation
            )
            push!(all_data, data)
        catch e
            println(e)
            println(string("Couldn't read ", file.name))
        end

    end

    sim_data = vcat(all_data...)    # join all dfs
    # filter(:variation => x -> x <= 15, all_data)
    
    return sim_data
end

# Read the human predictions
function read_subject_data(expt_id)

    folders = Dict(
        "Exp1" => "Exp1_allElasticities_fullMotion",
        "Exp2" => "Exp2_allElasticities_1second",
        "Exp3" => "Exp3_mediumElasticity_fullMotion",
        "Exp4" => "Exp4_mediumElasticity_1second"
    )
    exp_data_folder = joinpath(project_path, "Data", "HumanData", "EstimationTask", folders[expt_id], "Results")

    data = []
    for fname in readdir(exp_data_folder)
        sub_data = CSV.read(joinpath(exp_data_folder, fname), DataFrame)
        sub_data = insertcols(sub_data, :filename => fname)
        push!(data, sub_data) # this is a vector of dataframes
    end
    df = vcat(data...)     # combine all together with splat operator

    # elasticity is coded as integer
    id_no = parse(Int, expt_id[4])
    println("EXPT NO. ", id_no)
    if id_no >= 3
        df.elasticity = df.elasticity / 10
    end

    return df
end


################
# PROCESS DATA #
################

# Compute error on each stimulus, extract high error trials
function process_individual_stimuli_sim(sim_data)

    # Define model's elasticity judgment as average of output particles
    sim_data_pred = @chain sim_data begin
        @groupby :stimulusID
        @combine begin

            :judgment = mean(:elasticity)       
            :gtElasticity = first(:gtElasticity)  # :elasticity = :gtElasticity

            # compute error w.r.t. ground truth
            :error = mean(:elasticity) - first(:gtElasticity)
        end
        # @subset :elasticity .> 0.6
        @orderby :stimulusID
    end

    # Filter by error threshold, save high error trials
    high_error = sim_data_pred[sim_data_pred.error .> .2, :]
    high_error_trials = high_error[:, 1]

    return sim_data_pred, high_error_trials
end

# Group subject's elasticity judgments by stimulus and compute mean rating for each stimulus
function process_individual_stimuli_human(sub_data)

    nsubs = length(unique(sub_data.filename))   # number of subjects
    sub_data_pred = @chain sub_data begin
        @groupby :stimulusID
        # @groupby :elasticity# :filename
        # @DataFramesMeta.transform :prediction = mean(:rating)     # Gen also has a transform macro
        @combine begin

            :gtElasticity = mean(:elasticity)
            :judgment = mean(:rating)
            :std_err_mean = std(:rating) / sqrt(nsubs) 

        end
        # @subset :gtElasticity .> 0.6
        @orderby :stimulusID
    end

    return sub_data_pred
end


################
# DISPLAY DATA #
################

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


###########################
# GENERATING IDS AND PATH #
###########################

# Generate a string representation of the model's noise parameters, e.g. "Init00_Obs05_Tra05"
function generate_noise_id(args)

    init_string = split(string(args.init_vel_noise), ".")[2]
    if length(init_string) < 2
        init_string = string(init_string, "0")
    end

    obs_string = split(string(args.observation_noise), ".")[2]
    if length(obs_string) < 2
        obs_string = string(obs_string, "0")
    end

    transition_string = split(string(args.transition_noise), ".")[2] 
    if length(transition_string) < 2
        transition_string = string(transition_string, "0")
    end

    return string("Init", init_string, "_Obs", obs_string, "_Tra", transition_string)  
end

# Generate a string representation of the algorithm parameters
function generate_inference_param_id(args)

    particle_string = string(args.num_particles)
    rejuvenation_string = string(args.rejuvenation_moves)

    return string("Par", particle_string, "_Rej", rejuvenation_string)
end

# Generate the directory location where we will store the plots
function generate_plot_path(args, noise_id, inference_param_id)

    plots_path = joinpath(project_path, "Analysis")
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    plots_path = joinpath(plots_path, args.expt_id)
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    plots_path = joinpath(plots_path, args.model_id)
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    plots_path = joinpath(plots_path, args.target_id)
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    plots_path = joinpath(plots_path, noise_id)
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    plots_path = joinpath(plots_path, "SMC")
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    plots_path = joinpath(plots_path, inference_param_id)
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    plots_path = joinpath(plots_path, "Plots")
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    return plots_path
end
