# Utility functions to assist with file i/o for analysis scripts

using ZipFile
using CSV

project_path = dirname(@__DIR__)


###################
#  GENERATE PATH  #
###################

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

# Generate the directory location where we will store the plots
function generate_plot_path(expt_id, model_id, target_id, noise_id, type)

    plots_path = joinpath(project_path, "Analysis")
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    println("PLOTS PATH: ", plots_path)

    plots_path = joinpath(plots_path, expt_id)
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    println("PLOTS PATH: ", plots_path)

    plots_path = joinpath(plots_path, model_id)
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    plots_path = joinpath(plots_path, target_id)
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

    plots_path = joinpath(plots_path, "Plots")
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    plots_path = joinpath(plots_path, type)
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    return plots_path
end


##################
#  READING DATA  #
##################

# Read all simulation data files at the specified directory 
function read_simulation_data(expt_id, model_id, target_id, noise_id, algo, zip)
    
    all_data = []

    simulation_folder = joinpath(project_path, "Analysis", expt_id, model_id, target_id, noise_id, algo, "Data", "inferences")
    if zip 
        r = ZipFile.Reader(joinpath(simulation_folder, "inferences.zip")
        files = r.files
    else 
        files = readdir(simulation_folder)
    end

    for file in files
        
        # TODO: FINISH IMPLEMENTING NON-ZIP OPTIONs
        println(file.name)

        # Extract ground truth elasticity from filename
        elasticity_string, variation_string = split(file.name, ['_', '.'])
        elasticity = parse(Int, elasticity_string[end]) * 0.1 # fix order of magnitude
        variation = parse(Int64, variation_string[4:end])

        # Try to read the file into a dataframe
        try
            data = CSV.File(read(file)) |> DataFrame
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

    all_data = vcat(all_data...)    # join all dfs
    # filter(:variation => x -> x <= 15, all_data)
    
    return all_data
end


# Read the human predictions
function read_subject_data(expt_id)

    folders = Dict(
        1 => "Exp1_allElasticities_fullMotion",
        2 => "Exp2_allElasticities_1second",
        3 => "Exp3_mediumElasticity_fullMotion",
        4 => "Exp4_mediumElasticity_1second"
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
    if expt >= 3
        df.elasticity = df.elasticity / 10
    end

    return df
end
