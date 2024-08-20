# Utility functions to assist with fileio for analysis scripts

# Read the entire folder of simulation data
function read_simulation_data(expt, target_id)
    
    simulation_folder = joinpath(project_path, "Data", "BulletData", model_id, target_id, var_id, string("Exp", expt), "Inferences")

    all_data = []
    for file in readdir(simulation_folder)
        full_file_path = joinpath(simulation_folder, file)
        data = CSV.read(full_file_path, DataFrame)

        _, elasticity_string, variation = split(file, ['_', '.'])
        elasticity = parse(Int, elasticity_string[end]) * 0.1 # fix order of magnitude
    
        data = insertcols(
            data,
            "filename" => file,
            "stimulusID" => target_id * "_" * elasticity_string * "_" * variation,
            "gtElasticity" => elasticity,
            "variation" => parse(Int64, variation[4:end])
        )
        push!(all_data, data)
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