# Infer elasticity using MCMC

using Gen
using ZipFile

include("../Model/bouncing_object.jl")
include("../Utilities/fileio.jl")

PROPOSAL_STD_DEV = 0.5

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

# Resamples latent(s) from a Gaussian proposal distribution
@gen function proposal(trace::Gen.Trace)

    choices = get_choices(trace)

    # Resample restitution from a Gaussian centered at the previous estimate
    prev_res  = choices[:latents => 1 => :restitution]    
    println(prev_res)
    restitution = {:latents => 1 => :restitution} ~ trunc_norm(prev_res, 0.1, 0.0, 1.0)
    println(restitution)

    # prev_mass = choices[:latents => 1 => :mass]
    # mass = {:latents => 1 => :mass} ~ trunc_norm(prev_mass, 0.1, 0.0, Inf)

    return (restitution)
end

# Generate an initial trace then perform 1000 MH updates using Gaussian proposal
function gaussian_drift_inference(args, observations)

    # Generate an initial conditioned trace
    (trace, _) = generate(generate_trajectory, (args), observations)

    total_last_hundred = 0

    # Perform 1000 gaussian drift updates
    runs = 1000
    samples = 0
    accepts = 0 
    for i=1:runs        
        (trace, did_accept) = mh(trace, proposal, ())
        println(trace[:latents => 1 => :restitution], " Score: ", get_score(trace), " Did accept? ", did_accept)
        accepts += did_accept
        if i > runs - 100
            samples += 1 
            total_last_hundred += trace[:latents => 1 => :restitution]
        end
    end

    avg_last_hundred = total_last_hundred / 100
    println("Samples: ", samples)
    println("Accepts: ", accepts)
    
    return avg_last_hundred, trace
end