# Evaluate and plot model performance at inferring elasticity relative to ground truth
# Evaluate and plots model-inferred elasticities vs human judgments
#   Identify trials with poor model / human correlation by absolute error

using DataFrames
using DataFramesMeta
using Statistics
using Plots

include("../args.jl")
include("utils.jl")
include("plots.jl")


function main()

    args = Args()

    # Read simulation data into data frame
    noise_id = generate_noise_id(args)
    inference_param_id = generate_inference_param_id(args)
    sim_data = read_simulation_data(args.expt_id, args.model_id, args.target_id, noise_id, args.algorithm, inference_param_id, false)

    # Extract high-error trials
    sim_data, error_trials = process_individual_stimuli_sim(sim_data)
    println("High Error Trials: ", error_trials)


    # Generate plots

    # Generate output filepath for plots
    plots_path = generate_plot_path(args.expt_id, args.model_id, args.target_id, noise_id, inference_param_id, "TestJudgments")
    println("Plots Path: ", plots_path)

    # Set plot marker shape
    if args.target_id == "Cube"
        marker_shape = :square
    else
        marker_shape = :circle
    end

    # Plot model-inferred elasticity against ground truth
    plot_vs_gt("sim", sim_data, args.expt_id, args.target_id, marker_shape, plots_path)

    # If we have human data - RealFlow trials only
    if occursin("Exp", args.expt_id)

        # Read human data into data frame
        human_data = read_subject_data(args.expt_id)
        human_data = process_individual_stimuli_human(human_data)

        # Plot human judgments against ground truth
        plot_vs_gt("human", human_data, args.expt_id, args.target_id, marker_shape, plots_path)

        # Plot model judgments against human
        plot_sim_vs_human(sim_data, human_data, args.expt_id, args.target_id, marker_shape, plots_path)

    end

end

main()
