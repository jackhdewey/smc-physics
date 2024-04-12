# Inference and Prediction
#
# Defines a particle filter and proposal distribution that can be used to generate a set of concurrent simulations (particles)
# Each simulation is procedurally updated by one time step, scored, rejuvenated, and (potentially) resampled

include("../Model/bouncing_object.jl")
include("../Utilities/fileio.jl")


# Parses a sequence of observations from a choice map into a vector
function get_observations(choices::Gen.ChoiceMap, T::Int)
    observations = Vector{Gen.ChoiceMap}(undef, T)
    for i=1:T
        cm = choicemap()
        addr = :trajectory => i => :observation => :position
        set_submap!(cm, addr, get_submap(choices, addr))
        observations[i] = cm
    end
    return observations
end

# Resamples latent(s) from a Gaussian proposal distribution
@gen function proposal(trace::Gen.Trace)

    choices = get_choices(trace)

    # Resample restitution from a Gaussian centered at the previous estimate
    prev_res  = choices[:latents => 1 => :restitution]    
    restitution = {:latents => 1 => :restitution} ~ trunc_norm(prev_res, 0.1, 0.0, 1.0)

    # prev_mass = choices[:latents => 1 => :mass]
    # mass = {:latents => 1 => :mass} ~ trunc_norm(prev_mass, 0.1, 0.0, Inf)

    return (restitution)
end

# Particle filter
function infer(fname::String, id::String, gm_args::Tuple, obs::Vector{Gen.ChoiceMap}, num_particles::Int=20)

    # Extract trial identification    
    tokens = split(fname, "_")
    csv = split(tokens[3], ".")

    # Initiliaze the particle filter (with no observations)
    state = Gen.initialize_particle_filter(generate_trajectory, (gm_args[1], gm_args[2], 0), EmptyChoiceMap(), num_particles)
    
    # Iteratively simulate and filter particles
    argdiffs = (UnknownChange(), NoChange(), NoChange())
    for t=1:gm_args[3]
        @elapsed begin

            # Simulate the next time step and score against new observation
            Gen.particle_filter_step!(state, (gm_args[1],gm_args[2], t), argdiffs, obs[t])

            # Resample elasticity, resimulate to current time step, then accept / reject
            for i=1:num_particles
                state.traces[i], _ = Gen.mh(state.traces[i], proposal, ())
            end

            # Decide whether to cull and resample poorly performing particles
            Gen.maybe_resample!(state, ess_threshold=num_particles/2)
            
            # Dump current particles to a .csv file
            fname = string("BulletData/", id, "/Intermediate/particles_", tokens[2], "_", csv[1], "_", t, ".", csv[2])
            write_to_csv(state.traces, fname)

        end
    end

    # Gen.sample_unweighted_traces(state, 5)
    weights = Gen.get_log_weights(state)
    return state.traces, weights
end

# Given a set of trace (particles) with inferred elasticity, predict the next T postion ßobservations
function predict(particles, T::Int) 

    # For each particle
    n = length(particles)
    ppd = Vector{Gen.Trace}(undef, n)
    for i=1:n
        # Reset the number of time steps
        prev_args = get_args(particles[i])
        new_args = (prev_args[1], prev_args[2], prev_args[3] + T)
        arg_diffs = (UnknownChange(), NoChange(), NoChange())
        
        # Generate the next T time steps
        ppd[i] = first(Gen.update(particles[i], new_args, arg_diffs, choicemap()))
    end

    return ppd
end