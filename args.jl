# Arguments to the entire pipeline, separated out to ensure consistency / eliminate redundency between files

@with_kw struct Args

    # Data source
    expt_id::String = "BulletTest"
    gt_source::String = "Bullet"
    gt_shape::String = "Cube"       # Only used for Bullet

    # Model parameters
    model_id::String = "Modelv5"
    target_id::String = "Cube"
    noise_id::String = "PosVar075"
    init_vel_noise::Float32 = 0.0
    transition_noise::Float32 = 0.0
    observation_noise::Float32 = 0.0

    # Inference parameters
    algorithm::String = "SMC"       # MCMC, SMC, or DEBUG
    num_particles::Int = 20         # Only used for SMC
    save_intermediate::Bool = true  # Only used for SMC
    
    # Prediction parameters
    predict::Bool = false
    prediction_timesteps::Int = 90

    # Output filepath
    output_path::String = string(expt_id, "/", model_id, "/", target_id, "/", noise_id, "/", algorithm)

end