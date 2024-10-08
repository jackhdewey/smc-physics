# Helper functions for organizing, reading from, and writing data to .csv files

using DataFrames
using CSV

include("args.jl")

# Remove any filenames that contain certain tokens
function filter_unwanted_filenames(fnames)
    for i in ["predicted", "observed", "batch", "Store"]
        fnames = filter(!contains(i), fnames)
    end
    return fnames
end 

#######################
# READING AND WRITING #
#######################

# Extracts initial position, initial velocity, and trajectory from two .csv files
function read_obs_file(fname, pf::Bool=false)

    # Read ground truth initial state
    #fname = string("Data/RealFlowData/", fname) 
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

    # Split into a vector of choice maps 
    if pf
        observations = Vector{Gen.ChoiceMap}(undef, time_steps)

    # Or store the trajectory in a single choice map
    else 
        observations = Gen.choicemap()
    end

    for i=1:time_steps
        datum = values(trajectory_data[i, :])
        position = [datum[1], datum[2], datum[3]]
        addr = :trajectory => i => :observation => 1 => :position
        if pf
            observations[i] = Gen.choicemap((addr, position))
        else
            observations[addr] = position
        end
    end

    return initial_position, initial_orientation, initial_velocity, observations, time_steps
end

# Prepare directory for output files, and generate .zip file writers
function make_directories_and_writers(args::Args)
        
    # Locate / create base directory for output data
    dir_base = string("Data/BulletData/", args.output_id)
    if !isdir(dir_base)
        mkdir(dir_base)
    end

    # Locate / create directories for intermediate particles and output data
    particle_dir = string(dir_base, "/Intermediate/")
    if !isdir(particle_dir)
        mkdir(particle_dir)
    end

    output_dir = string(dir_base, "/Inferences/")
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    w1 = ZipFile.Writer(string(output_dir, "/inferences.zip"))
    w2 = ZipFile.Writer(string(particle_dir, "/particles.zip"))

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

############################################
# Sorting comparators (for ordering files) #
############################################

# Comparator to sort output files into correct order
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

# Comparator to sort intermediate particle filter state files into correct order
function trial_particle_order(x, y)

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

    as_int1 = parse(Int64, trial1)
    as_int2 = parse(Int64, trial2)

    if (as_int1 != as_int2)
        return as_int1 < as_int2
    end

    # Fourth token is particle number
    particle1 = replace(x_tokens[4], ".csv" => "")
    particle2 = replace(y_tokens[4], ".csv" => "")

    particle1_as_int = parse(Int64, particle1)
    particle2_as_int = parse(Int64, particle2)

    return particle1_as_int < particle2_as_int
end

# Comparator to sort plots of particle filter state over time into correct order
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

    # Third token is particle number
    particle1 = replace(x_tokens[3], ".png" => "")
    particle2 = replace(y_tokens[3], ".png" => "")

    particle1_as_int = parse(Int64, particle1)
    particle2_as_int = parse(Int64, particle2)

    return particle1_as_int < particle2_as_int
end
