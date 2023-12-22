# Provides functionality for visualizing the output of a particle filter

using Plots

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
    xlims --> (1, 120)
    yguide --> "height of cube (z)"
    ylims --> (0, 4.0)
    label --> false
    inds, z[inds, :]
end

function get_zs(trace::Gen.Trace)
    t, _... = get_args(trace)
    states = get_retval(trace)
    zs = Vector{Float64}(undef, t)
    for i = 1:t
        zs[i] = states[i].kinematics[1].position[3]
    end
    return zs
end

function animate_trace(trace::Gen.Trace; label="trace")
    t = first(get_args(trace))
    zs = reshape(get_zs(trace), (t, 1))
    @animate for i=2:t
        simplot(zs, i, label=label)
    end
end

function animate_traces(traces::Vector{<:Gen.Trace})
    zzs = reduce(hcat, map(get_zs, traces))
    t = size(zzs, 1)
    @animate for i=2:t
        simplot(zzs, i)
    end
end
