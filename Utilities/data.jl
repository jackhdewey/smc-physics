# Helper functions for reading from and writing data to a .csv file

# Extracts initial position, initial velocity, and trajectory from two .csv files
function read_observation_file(fname)

    # Read ground truth initial velocity
    fname = string("RealFlowData/", fname) 
    println("Reading...", fname)
    data = CSV.read(fname, DataFrame)
    datum = values(data[1, 11:13])
    initial_velocity = [datum[1], datum[2], datum[3]]

    # Read ground truth trajectory

    head, tail = split(fname, '.')
    fname2 = join([head, "_observed.", tail])
    println("Reading...", fname2)
    data = CSV.read(fname2, DataFrame)
    observations = Vector{Gen.ChoiceMap}(undef, size(data)[1])
    for i = 1:size(data)[1]
        addr = :trajectory => i => :observation => 1 => :position
        datum = values(data[i, :])
        new_datum = [datum[1], datum[3], datum[2]]
        cm = Gen.choicemap((addr, new_datum))
        observations[i] = cm
    end
    initial_position = get_value(observations[1], :trajectory => 1 => :observation => 1 => :position)

    return initial_position, initial_velocity, observations
end

function filter_unwanted_filenames(fnames)
    for i in ["predicted", "observed", "batch"]
        fnames = filter(!contains(i), fnames)
    end
    return fnames
end

function write_to_csv(particles, fname=joinpath(pwd(), "test.csv"))

    println("Writing simulation data to " * fname)
    particle_data = DataFrame(particle=Int[], elasticity=[], weight=[], frame=Int[], x=[], y=[], z=[])

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

    truncator(col, val) = trunc(val, digits=5)
    truncator(col, val::Int) = val

    CSV.write(fname, particle_data, transform=truncator)
end
