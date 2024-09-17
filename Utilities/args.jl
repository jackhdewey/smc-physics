@with_kw struct Args

    # Data source
    gt_source::String = "Bullet"

    # Model parameters
    model_id::String = "Modelv5"
    target_id::String = "Cube"
    noise_id::String = "PosVar05"
    expt_id::String = "BulletxBullet"
    output_id::String = string(model_id, "/", target_id, "/", noise_id, "/", expt_id)

    # Inference parameters
    algorithm::String = "MCMC"    # MCMC, PARTICLE_FILTER, or DEBUG
    num_particles::Int = 20
    save_intermediate::Bool = true

    # Prediction parameters
    predict::Bool = false
    prediction_timesteps::Int = 90

end
