# Arguments to the entire pipeline, separated out to ensure consistency / eliminate redundency between files

@with_kw struct Args

    # Data source
    gt_source::String = "Bullet"
    gt_shape::String = "Cube"       # Only used for PyBullet
    expt_id::String = "BulletTest"

    # Model parameters
    model_id::String = "Modelv5"
    target_id::String = "Cube"
    noise_id::String = "PosVar075"

    # Inference parameters
    algorithm::String = "SMC"       # MCMC, SMC, or DEBUG
    num_particles::Int = 20
    save_intermediate::Bool = true
    
    # Prediction parameters
    predict::Bool = false
    prediction_timesteps::Int = 90

    # Output filepath
    output_path::String = string(model_id, "/", target_id, "/", noise_id, "/", expt_id, "/", algorithm)

end