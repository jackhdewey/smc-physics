#########
# PLOTS #
#########

# Plot model against ground truth
function plot_vs_gt(type, data, expt_id, target_id, marker_shape, plots_path)

    if type == "sim"
        ylabel = "Model Estimate"
        title = "$expt_id Model ($target_id) v. Ground Truth Individual Stimuli"
    else 
        ylabel = "Human Estimate"
        title = "$expt_id Human v. Ground Truth Individual Stimuli"
    end

    p = palette(:jet)
    scatter(data.gtElasticity,
            data.judgment,
            #yerror=sim.std_err_mean ./ 2,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor = data.gtElasticity,  
            clims=(0, 1),
            xlims=(0, 1.05),
            ylims=(0, 1),
            xlabel = "Ground Truth Elasticity",
            ylabel = ylabel,
            title = title,
            colorbar=true,
            legend=false,
            palette=p,
            markershape = marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(data.gtElasticity, data.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)

    type == "sim" ?
        savefig(joinpath(plots_path, "sim_vs_gt_individual_stimuli.png")) :
        savefig(joinpath(plots_path, "human_vs_gt_individual_stimuli.png"))

end


# Plot mean model estimates against mean human judgments
function plot_sim_vs_human(sim_data, human_data, expt_id, target_id, marker_shape, plots_path)

    yerror = human_data.std_err_mean ./ 2

    title = "$expt_id Model ($target_id) vs. Human Individual Stimuli"
    filename = "sim_vs_human_individual_stimuli"

    # @infiltrate
    # default(aspect_ratio = :equal)
    p = palette(:jet)
    scatter(human_data.judgment,
            sim_data.judgment,
            yerror = yerror,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor = human_data.gtElasticity,
            xlims=(0, 1.05),
            ylims=(0, 1),
            clims=(0, 1),
            xlabel="Human",
            ylabel="Model",
            title=title,
            legend=false,
            colorbar=true,
            palette=p,
            markershape=marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(sim_data.judgment, human_data.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)

    savefig(joinpath(plots_path, filename))


    title = "$expt_id Model ($target_id) vs. Human Mean"
    filename = "sim_vs_human_average_judgments"

    sim_data = @chain sim_data begin
        @groupby :gtElasticity
        @combine :judgment = mean(:judgment)
    end

    human_data = @chain human_data begin
        @groupby :gtElasticity
        @combine :judgment = mean(:judgment)
    end
    # gui()

    # @infiltrate
    # default(aspect_ratio = :equal)
    p = palette(:jet)
    scatter(human_data.judgment,
            sim_data.judgment,
            yerror = yerror,
            aspect_ratio=:equal,
            markersize=5,
            markeralpha=0.5,
            zcolor = human_data.gtElasticity,
            xlims=(0, 1.05),
            ylims=(0, 1),
            clims=(0, 1),
            xlabel="Human",
            ylabel="Model",
            title=title,
            legend=false,
            colorbar=true,
            palette=p,
            markershape=marker_shape)

    plot!(0:1, 0:1, line=:dash)
    corr_string = "r = " * string(round(cor(sim_data.judgment, human_data.judgment), digits=3))
    annotate!(0.2, 0.8, corr_string, 10)

    savefig(joinpath(plots_path, filename))

end

# Plot all particle trajectories 
function plot_trajectories(gt_dir, gt_files, data_path, particle_files, particle_indices, num_particles, plot_path)

    # Plot parameters
    if interactive
        pyplot()
    end

    for i in eachindex(gt_files)

        # Read ground truth file to dataframe
        gt_file = string(gt_dir, gt_files[i])
        ground_truth = CSV.read(gt_file, DataFrame)

        # Generate plot base
        tokens = split(gt_files[i], "_")
        title = string("Stimulus: ", tokens[2], " ", tokens[3], "\n")
        plt = plot3d(
                1,
                xlim=(-0.5, 0.5),
                ylim=(-0.5, 0.5),
                zlim=(0, 1),
                title=title,
                legend=false,
                marker=2,
                seriestype=:scatter,
                size=(1200, 800),
                gridlinewidth=8
        )

        # Procedurally generate plot
        true_x = []
        true_y = []
        true_z = []
        num_timesteps = size(ground_truth)[1]
        for t = 1:num_timesteps

            # Extend ground truth trajectory by one time step and add to plot
            true_x = [true_x; ground_truth[t, 1]]
            true_y = [true_y; ground_truth[t, 2]]
            true_z = [true_z; ground_truth[t, 3]]
            plot!(plt, true_x, true_y, true_z, linewidth=3, linecolor=:red)

            # Extract and sort particle files
            #files = map((file) -> file.name, reader.files)

            # Index into correct particle file and read into dataframe
            d_file = string(data_path, "/", particle_files[particle_indices[i] + t - 1])
            #t_file = filter((file) -> file.name == d_file, reader.files)[1]
            data = CSV.File(read(d_file)) |> DataFrame

            # For each particle
            for p = 1:num_particles
                particle = data[data.particle.==p, :]

                # elasticity = data[data.particle .== i, 2]
                # weight = data[data.particle .== i, 3]
                # @df plot!(plt1, particle[:, 5:7])

                # Read particle trajectory up to current time step
                x_trajectory = []
                y_trajectory = []
                z_trajectory = []
                for row = 1:t
                    x_trajectory = [x_trajectory; particle[row, 5]]
                    y_trajectory = [y_trajectory; particle[row, 6]]
                    z_trajectory = [z_trajectory; particle[row, 7]]
                end

                # Add particle trajectory to plot
                # TODO: Set line color / alpha based on elasicity 
                plot!(plt, x_trajectory, y_trajectory, z_trajectory)

            end

            # Every fifth time step, display the plot and (if static) save as .png
            if t % plot_interval == 0

                # Update title
                title!(string(title, "\nTimestep: ", t, " / ", num_timesteps))

                display(plt)

                directories = split(d_file, "/")
                println(directories[10])
                tokens = split(directories[10], "_")
                if !interactive
                    savefig(string(plot_path, "/", tokens[1], "_", tokens[2], "_", t))
                end
            end
        end
    end
end
