# Arguments to the entire pipeline, separated out to ensure consistency + eliminate redundency between files

using Parameters

@with_kw mutable struct Args

    # Data source
    expt_id::String = "BulletTest"
    gt_source::String = "Bullet"     # Bullet or RealFlow
    gt_shape::String = "Cube"        # Only used for Bullet

    # Model parameters
    model_id::String = "Modelv5"
    target_id::String = "Cube"
    init_vel_noise::Float32 = 0.0
    observation_noise::Float32 = 0.05
    transition_noise::Float32 = 0.075

    # Inference parameters
    algorithm::String = "SMC"        # MCMC, SMC, or DEBUG
    num_particles::Int = 20          # Only used for SMC
    rejuvenation_moves::Int = 5      # Only used for SMC
    save_particles::Bool = false     # Only used for SMC
    
    # Prediction parameters
    predict::Bool = false
    prediction_timesteps::Int = 90

    # Analysis parameters
    plot_mean = true

end