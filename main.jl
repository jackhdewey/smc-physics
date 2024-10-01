# Galileo 3
#
# Infers the elasticity of a bouncing soft body from a sequence of postition observations and predicts its future trajectory 
# What representation in the generative model produces inferences that correlate best with human judgments?
#      * using a rigid body
#      * using a sphere
#      * etc.
# Ground truth trajectories are generated in RealFlow and read from .csv files
#
# TODO: Infer elasticity for trajectories simulated in PyBullet
# TODO: Debug MCMC
# TODO: Display particle trajectory data 

using Distributed
using Parameters
using ZipFile

include("Model/bouncing_object.jl")
include("Inference/mcmc.jl")
include("Inference/particle_filter.jl")
include("Utilities/plots.jl")


@with_kw struct Args

    # Data source
    gt_source::String = "Bullet"

    # Experiment
    expt_id::String = "BulletTest"

    # Model parameters
    model_id::String = "Modelv5"
    target_id::String = "Sphere"
    noise_id::String = "PosVar01"

    # Filepaths
    output_id::String = string(model_id, "/", target_id, "/", noise_id, "/", expt_id)

    # Inference parameters
    algorithm::String = "MCMC"    # MCMC, PARTICLE_FILTER, or DEBUG
    num_particles::Int = 20
    save_intermediate::Bool = true

    # Prediction parameters
    predict::Bool = false
    prediction_timesteps::Int = 90

end

@everywhere begin
function run(fname, args, w1, w2)
    
    # Initialize simulation context 
    if args.algorithm=="DEBUG"
        client = bullet.connect(bullet.GUI)::Int64
    else 
        client = bullet.connect(bullet.DIRECT)::Int64
    end
    
    bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
    bullet.resetDebugVisualizerCamera(4, 0, -4, [0, 0, 2])

    if gt_source=="Bullet"
        sim = BulletSim(step_dur=1/60; client=client)
    else 
        sim = BulletSim(step_dur=1/30; client=client)
    end

    # Initialize scene 
    init_scene(args.algorithm=="DEBUG")

    # Initialize target object state using observation data
    args.gt_source == "RealFlow" ? 
        fname = string("Data/RealFlowData/", args.expt_id, "/", fname) : 
        fname = string("Tests/BulletStimulus/Data/", args.target_id, "/", fname)
    init_position, _, init_velocity, observations, t_s = read_obs_file(fname, args.algorithm=="PARTICLE_FILTER")
    init_state = init_target_state(sim, args.target_id, init_position, init_velocity)
    model_args = (sim, init_state, t_s)

    if args.algorithm == "DEBUG"

        # Generate and display a single trace 
        println("DEBUGGING")
        (trace, _) = Gen.generate(generate_trajectory, model_args, observations)
        display(Gen.get_choices(trace))

    elseif args.algorithm == "MCMC"
        
        println("Initializing MCMC")
        avg_last_hundred, trace = gaussian_drift_inference(model_args, observations)
        #println(trace[:latents => 1 => :restitution])
        println("Elasticity estimate: ", avg_last_hundred)

        # Write output trajectory to a .csv file
        tokens = split(fname, "_")
        tokens = split(fname, "_")
        fname = string(tokens[2], "_", tokens[3])
        println(fname)
        #f = ZipFile.addfile(w1, fname)
        #write_to_csv(results, f)
        
    else

        # Filter n particles to explain the complete trajectory
        println("Initializing Particle Filter")
        results, _ = infer(generate_trajectory, model_args, observations, w2, args.num_particles, args.save_intermediate, fname)

        # Write output particles to a .csv file
        tokens = split(fname, "_")
        fname = string(tokens[2], "_", tokens[3])
        println(fname)
        f = ZipFile.addfile(w1, fname)
        write_to_csv(results, f)

        if args.predict

            # For each particle, simulate the next 90 time steps
            ppd = predict(results, args.prediction_timesteps)
            #gif(animate_traces(ppd), fps=24)
        
            # Write predicted trajectories to a .csv file
            fname = string(dir_base, "/Predictions/predictions_", tokens[2], "_", tokens[3])
            write_to_csv(ppd, fname)

        end
    end

    bullet.disconnect()
    #=
    catch e
        println(e) 
        println("Disconnecting Bullet")
        bullet.disconnect()
    end
    =#
end
end

########
# MAIN #
########

function main()

    args = Args()

    # Read ground truth trajectories
    args.gt_source == "RealFlow" ? 
        dir = string("Data/RealFlowData/", args.expt_id, "/") : 
        dir = string("Tests/BulletStimulus/Data/", args.target_id, "/")
    fnames = readdir(dir)
    fnames = filter_unwanted_filenames(fnames)
    sort!(fnames, lt=trial_order)

    w1, w2 = make_directories_and_writers(args.output_id)

    run(fnames[1], args, w1, w2)
    #=
    for fname in fnames
        run(fname, args, w1, w2)
    end
    =#

    # Map each observed trajectory to a process executing a particle filter
    #pmap(fname -> run(fname, args, w1, w2), fnames)

    close(w1)
    close(w2)  
end

main()
