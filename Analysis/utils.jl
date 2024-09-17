# Utility functions to assist with file i/o for analysis scripts

using ZipFile
using CSV


# Read the entire folder of simulation data
function read_simulation_data(model_id, target_id, noise_id, expt_id)
    
    simulation_folder = joinpath(project_path, "Data", "BulletData", model_id, target_id, noise_id, expt_id, "Inferences")
    r = ZipFile.Reader(joinpath(simulation_folder, "inferences.zip"))

    all_data = []

    #=
    for file in readdir(simulation_folder)
        full_file_path = joinpath(simulation_folder, file)
        data = CSV.read(full_file_path, DataFrame)

        _, elasticity_string, variation = split(file, ['_', '.'])
        elasticity = parse(Int, elasticity_string[end]) * 0.1 # fix order of magnitude
    =#

    for file in r.files
        println(file.name)

        # Extract ground truth elasticity from filename
        elasticity_string, variation_string = split(file.name, ['_', '.'])
        elasticity = parse(Int, elasticity_string[end]) * 0.1 # fix order of magnitude
        variation = parse(Int64, variation_string[4:end])

        #=
        if file.name == "Ela4_Var1.csv" || file.name == "Ela7_Var1.csv"
            continue
        end
        =#
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
        catch
            println(string("Couldn't read", file.name))
        end

    end

    all_data = vcat(all_data...)    # join all dfs
    # filter(:variation => x -> x <= 15, all_data)
    
    return all_data
end


# Read the human predictions
function read_subject_data(expt)

    folders = Dict(
        1 => "Exp1_allElasticities_fullMotion",
        2 => "Exp2_allElasticities_1second",
        3 => "Exp3_mediumElasticity_fullMotion",
        4 => "Exp4_mediumElasticity_1second"
    )
    exp_data_folder = joinpath(project_path, "Data", "HumanData", "EstimationTask", folders[expt], "Results")

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

# Generate the filepath where we will save the plots
function generate_plot_path(model_id, target_id, noise_id, expt_id)

    plots_path = joinpath(project_path, "Analysis", "Plots")
    if !isdir(plots_path)
        mkdir(plots_path)
    end

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

    plots_path = joinpath(plots_path, expt_id)
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    plots_path = joinpath(plots_path, "Judgments")
    if !isdir(plots_path)
        mkdir(plots_path)
    end

    return plots_path
end