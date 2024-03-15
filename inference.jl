# Inference and Prediction

include("elasticity.jl")
include("Utilities/fileio.jl")

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
function filter(gm_args::Tuple, obs::Vector{Gen.ChoiceMap}, num_particles::Int=20, fname)

    # Extract trial identification    
    tokens = split(fname, "_")
    csv = split(tokens[3], ".")

    # Initiliaze the particle filter (with no observations)
    state = Gen.initialize_particle_filter(generate_trajectory, (0, gm_args[2:3]...), EmptyChoiceMap(), num_particles)
    
    # Iteratively simulate and filter particles
    argdiffs = (UnknownChange(), NoChange(), NoChange())
    for (t, obs) = enumerate(obs)
        @elapsed begin
            
            # Decide whether to cull and resample poorly performing particles
            Gen.maybe_resample!(state, ess_threshold=num_particles/2)

            # Simulate the next time step and score against new observation
            Gen.particle_filter_step!(state, (t, gm_args[2:3]...), argdiffs, obs)

            # Resample elasticity, resimulate to current time step, then accept / reject
            for i=1:num_particles
                state.traces[i], _ = Gen.mh(state.traces[i], proposal, ())
            end
            
            # Dump current particles to a .csv file
            fname = string("BulletData/Intermediate/particles_", tokens[2], "_", csv[1], "_", t, ".", csv[2])
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
        new_args = (prev_args[1] + T, prev_args[2:end]...)
        arg_diffs = (UnknownChange(), NoChange(), NoChange())
        
        # Generate the next T time steps
        ppd[i] = first(Gen.update(particles[i], new_args, arg_diffs, choicemap()))
    end

    return ppd
end