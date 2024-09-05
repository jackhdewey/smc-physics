# Tests whether a range of elasticity settings produce plausible trajectories
function test_elasticity(sim, init_state)
    
    gt_constraints = choicemap((:latents => 1 => :restitution, 0.2))

    args = (60, init_state, sim)
    trace = first(generate(generate_trajectory, args, gt_constraints))

    traces = [trace]
    for i=2:5
        res = i * 0.2
        gt_constraints = choicemap((:latents => 1 => :restitution, res))
        trace = first(generate(generate_trajectory, args, gt_constraints))
        push!(traces, trace)
    end

    gif(animate_traces(traces), fps=24)

    for trace in traces
        choices = get_choices(trace)
        positions = [choices[:trajectory => i => :observation => 1 => :position] for i=10:20]
        ys = map(pos -> pos[2], positions)
        zs = map(pos -> pos[3], positions)
        y_max = maximum(ys)
        z_max = maximum(zs)
        println("y_max: ", y_max)
        println("z_max: ", z_max)
        coefficient_of_restitution = z_max / 0.9
    end
end