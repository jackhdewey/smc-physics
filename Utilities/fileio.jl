# Helper functions for organizing, reading from, and writing data to .csv files

using DataFrames
using CSV


#############################
# READING OBSERVATION FILES #
#############################

# Removes unwanted filenames from a directory by checking for tokens
function filter_unwanted_filenames(fnames)
    for i in ["predicted", "observed", "batch", "Store"]
        fnames = filter(!contains(i), fnames)
    end
    return fnames
end 

# Extracts initial position, initial velocity, and trajectory from two .csv files
function read_obs_file(fname, alg::String)

    # Read ground truth initial state
    println("Reading...", fname)
    init_state = CSV.read(fname, DataFrame)

    # Initial position
    datum = values(init_state[1, 5:7])
    initial_position = [datum...]

    # Initial orientation
    datum = values(init_state[1, 8:10])
    initial_orientation = [datum...]

    # Initial velocity
    datum = values(init_state[1, 11:13])
    initial_velocity = [datum...]

    # Read ground truth trajectory
    head, tail = split(fname, '.')
    obs_fname = join([head, "_observed.", tail])
    println("Reading...", obs_fname)
    trajectory_data = CSV.read(obs_fname, DataFrame)
    time_steps = size(trajectory_data)[1]

    if alg == "SMC"     # If SMC, split into a vector of choice maps 
        observations = Vector{Gen.ChoiceMap}(undef, time_steps)
    else                # Otherwise store the full trajectory in a single choice map
        observations = Gen.choicemap()
    end

    # Populate choice map(s)
    for i=1:time_steps
        datum = values(trajectory_data[i, :])
        position = [datum[1], datum[2], datum[3]]
        addr = :trajectory => i => :observation => 1 => :position
        if alg == "SMC"
            observations[i] = Gen.choicemap((addr, position))
        else
            observations[addr] = position
        end
    end

    return initial_position, initial_orientation, initial_velocity, observations, time_steps
end


################################
# Comparators for File Sorting #
################################

# Sort output files based on elasticity -> variation
function trial_order(x, y)

    x_tokens = split(x, "_")
    y_tokens = split(y, "_")

    # Second token is elasticity
    elasticity1 = replace(x_tokens[2], "Ela" => "")
    elasticity2 = replace(y_tokens[2], "Ela" => "")

    as_int1 = parse(Int64, elasticity1)
    as_int2 = parse(Int64, elasticity2)

    if (as_int1 != as_int2)
        return as_int1 < as_int2
    end

    # Third token is trial number
    trial1 = replace(x_tokens[3], "Var" => "")
    trial2 = replace(y_tokens[3], "Var" => "")

    trial1 = replace(trial1, ".csv" => "")
    trial2 = replace(trial2, ".csv" => "")

    as_int1 = parse(Int64, trial1)
    as_int2 = parse(Int64, trial2)

    if (as_int1 != as_int2)
        return as_int1 < as_int2
    end
end

# Sort intermediate particle filter state files baseed on elasticity -> variation -> time step
function trial_particle_order(x, y)

    x_tokens = split(x, "_")
    y_tokens = split(y, "_")

    # First token is elasticity
    elasticity1 = replace(x_tokens[1], "Ela" => "")
    elasticity2 = replace(y_tokens[1], "Ela" => "")

    as_int1 = parse(Int64, elasticity1)
    as_int2 = parse(Int64, elasticity2)

    if (as_int1 != as_int2)
        return as_int1 < as_int2
    end

    # Second token is trial number
    trial1 = replace(x_tokens[2], "Var" => "")
    trial2 = replace(y_tokens[2], "Var" => "")

    as_int1 = parse(Int64, trial1)
    as_int2 = parse(Int64, trial2)

    if (as_int1 != as_int2)
        return as_int1 < as_int2
    end

    # Third token is time step
    step1 = replace(x_tokens[3], ".csv" => "")
    step2 = replace(y_tokens[3], ".csv" => "")

    as_int1 = parse(Int64, step1)
    as_int2 = parse(Int64, step2)

    return as_int1 < as_int2
end

# Sort plots of particle filter state based on elasticity -> variation -> time step
function png_particle_order(x, y)

    x_tokens = split(x, "_")
    y_tokens = split(y, "_")

    # First token is elasticity
    elasticity1 = replace(x_tokens[1], "Ela" => "")
    elasticity2 = replace(y_tokens[1], "Ela" => "")

    as_int1 = parse(Int64, elasticity1)
    as_int2 = parse(Int64, elasticity2)

    if (as_int1 != as_int2)
        return as_int1 < as_int2
    end

    # Second token is trial number
    trial1 = replace(x_tokens[2], "Var" => "")
    trial2 = replace(y_tokens[2], "Var" => "")

    as_int1 = parse(Int64, trial1)
    as_int2 = parse(Int64, trial2)

    if (as_int1 != as_int2)
        return as_int1 < as_int2
    end

    # Third token is time step
    step1 = replace(x_tokens[3], ".png" => "")
    step2 = replace(y_tokens[3], ".png" => "")

    as_int1 = parse(Int64, step1)
    as_int2 = parse(Int64, step2)

    return as_int1 < as_int2
end


########################
# WRITING OUTPUT FILES #
########################

# Generate a string representation of the noise parameters for a model variation, e.g. "Init00_Obs05_Tra05"
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

# Prepare directory for output files, and generate .zip file writers
function make_directories(output_path)
    
    tokens = split(output_path, "/")

    # If it doesn't already exist, create base directory for output data
    expt_dir = string("Analysis/", tokens[1], "/")
    if !isdir(expt_dir)
        mkdir(expt_dir)
    end
    model_dir = string(expt_dir, tokens[2], "/")
    if !isdir(model_dir)
        mkdir(model_dir)
    end
    target_dir = string(model_dir, tokens[3], "/")
    if !isdir(target_dir)
        mkdir(target_dir)
    end
    noise_dir = string(target_dir, tokens[4], "/")
    if !isdir(noise_dir)
        mkdir(noise_dir)
    end
    alg_dir = string(noise_dir, tokens[5], "/")
    if !isdir(alg_dir)
        mkdir(alg_dir)
    end
    infer_dir = string(alg_dir, tokens[6], "/")
    if !isdir(infer_dir)
        mkdir(infer_dir)
    end
    data_dir = string(infer_dir, "Data/")
    if !isdir(data_dir)
        mkdir(data_dir)
    end

    return data_dir
end

function make_writers(data_dir, algorithm)

    # Initialize zip writers for inference data and (if SMC) intermediate particles
    w1 = ZipFile.Writer(string(data_dir, "/inferences.zip"))
    w2 = nothing
    if algorithm == "SMC"
        w2 = ZipFile.Writer(string(data_dir, "/particles.zip"))
    end
    
    return w1, w2
end

# Writes data for each particle (elasticity, log weight, and trajectory) to a .csv file
function write_to_csv(particles, fname=joinpath(pwd(), "test.csv"))

    # Initialize a data frame with appropriate 
    #println("Writing simulation data to " * fname)
    particle_data = DataFrame(particle=Int[], elasticity=[], weight=[], frame=Int[], x=[], y=[], z=[])

    # Iterate over the particles 
    for (p, particle) in enumerate(particles)
        ela = particle[:latents => 1 => :restitution]
        for (f, frame) in enumerate(particle[:trajectory])
            pos = convert(Vector, frame.kinematics[1].position)
            #ori = convert(Vector, frame.kinematics[1].orientation)
            weight = get_score(particle)
            data = [p; ela; weight; f; pos]
            push!(particle_data, data)
        end
    end

    # Truncate to 5 digits
    truncator(col, val) = trunc(val, digits=5)
    truncator(col, val::Int) = val

    CSV.write(fname, particle_data, transform=truncator)
end

function zip_folder(filepath)

    w1, w2 = make_writers(filepath)

    inference_path = string(filepath, "/inferences/")
    particle_path = string(filepath, "/intermediate/")
    println(inference_path)
    println(particle_path)

    for (root, dirs, files) in walkdir(filepath)
    end

end
