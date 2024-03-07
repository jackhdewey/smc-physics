# Provides functionality for visualizing the output of a particle filter

using Plots
using Gen

# Generates an animated 2D plot of height over time 
@userplot SimPlot
@recipe function f(cp::SimPlot)
    z, t = cp.args
    cs = size(z, 1)
    k = 10
    inds = (max(1, t-k):t)
    n = length(inds)
    linewidth -->range(0,10, length=n)
    seriesalpha --> range(0,1,length=n)
    xguide --> "time"
    xlims --> (1, 40)
    yguide --> "height of cube (z)"
    ylims --> (0, 1.0)
    label --> false
    inds, z[inds, :]
end

# Extract height measurements from the given trace
function get_zs(trace::Gen.Trace)
    t, _... = get_args(trace)
    states = get_retval(trace)
    zs = Vector{Float64}(undef, t)
    for i = 1:t
        zs[i] = states[i].kinematics[1].position[3]
    end
    return zs
end

# Animate a single trajectory
function animate_trace(trace::Gen.Trace; label="trace")
    t = first(get_args(trace))
    zs = reshape(get_zs(trace), (t, 1))
    @animate for i=2:t
        simplot(zs, i, label=label)
    end
end

# Animate multiple trajectories in one plot
function animate_traces(traces::Vector{<:Gen.Trace})
    zzs = reduce(hcat, map(get_zs, traces))
    t = size(zzs, 1)
    @animate for i=2:t
        simplot(zzs, i)
    end
end

# Plots

# Plot a single execution trace (particle) of the generative model
function visualize_particle_unfold(trace; show_data=true, max_T=get_args(trace)[1], overlay=false)

    # Extract the number of time steps and complete state sequence
    (T,) = Gen.get_args(trace)
    choices = Gen.get_choices(trace)
    (init_state, states) = Gen.get_retval(trace)

    # Populate vectors with observations and estimated positions
    xs = Vector{Float64}(undef, T+1)
    ys = Vector{Float64}(undef, T+1)
    obs = Vector{Float64}(undef, T+1)
    obs[1] = choices[:init_state => :obs => :bearing]
    xs[1] = init_state.x
    ys[1] = init_state.y
    for i=1:T
        obs[i+1] = choices[:trajectory => i => :obs => :bearing]
        xs[i+1] = states[i].x
        ys[i+1] = states[i].y
    end

    # Plot the estimated positions 
    f = overlay ? scatter! : scatter
    fig = f(xs[1:max_T+1], ys[1:max_T+1], s=:auto, label=nothing)

    # Plot the estimated bearings
    if show_data
        for z in obs[1:max_T+1]
            dx = cos(z) * 0.5
            dy = sin(z) * 0.5
            plot!([0., dx], [0., dy], color="red", alpha=0.3, label=nothing)
        end
    end
end

# Overlay multiple particles onto a single plot
function overlay_particles(renderer, traces; same_data=true, args...)

    fig = plot(title="Observed bearings (red) and \npositions of individual traces (one color per trace)", xlabel="X", ylabel="Y")
    
    renderer(traces[1], show_data=true, overlay=true, args...)
    for i=2:length(traces)
        renderer(traces[i], show_data=!same_data, overlay=true, args...)
    end

    return fig
end

# 
function plot_traj(all_particles)
    traj = filter(:particle => x -> x == 1, all_particles)
    plt = plot3d(traj.x, traj.y, traj.z)
end
