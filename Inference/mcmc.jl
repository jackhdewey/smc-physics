# Infer elasticity using MCMC

using Gen
using ZipFile

include("../Model/bouncing_object.jl")
include("../Utilities/fileio.jl")


# Propose a new elasticity by resampling, resimulating, and MH accept/reject
function block_resimulation_update(trace)
    
    # Select elasticity, then resample, test, and accept/reject
    params = Gen.select(:latents => 1 => :restitution)
    (trace, _) = metropolis_hastings(trace, params)

    return trace
end

# Generate an initial trace then perform 500 block resimulation updates
function block_resimulation_inference(args, observations)

    # Generate an initial conditioned trace
    (trace, _) = generate(generate_trajectory, (args), observations)

    # Perform 500 block resimulation updates
    for _=1:500
        trace = block_resimulation_update(trace)
    end

    return trace
end

# Resamples latent(s) from a Gaussian proposal distribution - for use in MCMC rejuvenation
@gen function proposal(trace::Gen.Trace)

    choices = get_choices(trace)

    # Resample restitution from a Gaussian centered at the previous estimate
    prev_res  = choices[:latents => 1 => :restitution]    
    restitution = {:latents => 1 => :restitution} ~ trunc_norm(prev_res, 0.1, 0.0, 1.0)

    # prev_mass = choices[:latents => 1 => :mass]
    # mass = {:latents => 1 => :mass} ~ trunc_norm(prev_mass, 0.1, 0.0, Inf)

    return (restitution)
end

# Generate an initial trace then perform 500 block resimulation updates
function gaussian_drift_inference(args, observations)

    # Generate an initial conditioned trace
    (trace, _) = generate(generate_trajectory, (args), observations)

    # Perform 500 block resimulation updates
    for _=1:1000
        (trace, did_accept) = mh(trace, proposal, ())
    end

    return trace
end