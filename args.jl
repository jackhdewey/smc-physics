# Arguments to the entire pipeline, separated out to ensure consistency + eliminate redundency between files

using Parameters

@with_kw mutable struct Args

    # Data source
    expt_id::String = "Exp1"
    gt_source::String = "RealFlow"
    gt_shape::String = "Cube"       # Only used for Bullet

    # Model parameters
    model_id::String = "Modelv5"
    target_id::String = "Sphere"
    init_vel_noise::Float32 = 0.0
    observation_noise::Float32 = 0.05
    transition_noise::Float32 = 0.075

    # Inference parameters
    algorithm::String = "SMC"       # MCMC, SMC, or DEBUG
    num_particles::Int = 20         # Only used for SMC
    save_particles::Bool = false    # Only used for SMC
    
    # Prediction parameters
    predict::Bool = false
    prediction_timesteps::Int = 90

end