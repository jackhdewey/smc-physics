module PhyBullet

using PyCall
using Parameters
using StaticArrays
using DocStringExtensions

export pb,
    BulletSim,
    BulletState,
    BulletElement,
    BulletElemState,
    BulletElemLatents

const pb = PyNULL()
function __init__()
    copy!(pb, pyimport("pybullet"))
end

"""
Pass along abstract types to be implemented by e.g. rigid_body.jl  
"""
abstract type BulletElement <: Element{BulletSim} end
abstract type BulletElemState{T<:BulletElement} <: ElemState{T} end
abstract type BulletElemLatents{T<:BulletElement} <: ElemLatents{T} end

function get_state end

function set_state! end

function get_latents end

function set_latents! end

include("Elements/rigid_body.jl")


"""
Global simulation context for the Bullet physics engine
"""
@with_kw struct BulletSim <: PhySim
    # Client id for pybullet
    client::Int64
    # Timestep duration of bullet engine (default: 4.2ms)
    pb_timestep::Float64 = 1 / 240
    # Amount of time between `forward_steps` (default=16.7ms)
    step_dur::Float64 = 1 / 60
end

"""
State of a Bullet simulation
TODO: Remove latents and kinematics, add collisions
"""
struct BulletState <: PhyState{BulletSim}
    elements::AbstractVector{BulletElement}
    latents::AbstractVector{BulletElemLatents}
    kinematics::AbstractVector{BulletElemState}
end

# Constructs a bullet state with the provided elements
function BulletState(sim::BulletSim, elements::AbstractVector{T}) where {T<:BulletElement}
    latents = map(x -> get_latents(x, sim), elements)
    kinematics = map(x -> get_state(x, sim), elements)
    BulletState(elements, latents, kinematics)
end

function PhySMC.sync!(sim::BulletSim, world_state::BulletState)
    for (elem, ls, est) in zip(world_state.elements,
                               world_state.latents,
                               world_state.kinematics)
        set_state!(elem, sim, est)
        set_latents!(elem, sim, ls)
    end
    return nothing
end

# Run the physics engine forward by one time step to generate the next state
function PhySMC.forward_step(sim::BulletSim, state::BulletState)
    
    # progress by `st.step_dur`
    dt::Float64 = 0.0
    while dt <= sim.step_dur
        @pycall pb.stepSimulation(;
                                  physicsClientId = sim.client
                                  )::PyObject
        dt += sim.pb_timestep
    end

    # extract resulting state
    ne = length(st.elements)
    elements = Vector{BulletElement}(undef, ne)
    latents = Vector{BulletElemLatents}(undef, ne)
    kinematics = Vector{BulletElemState}(undef, ne)
    @inbounds for i = 1:ne
        elements[i] = state.elements[i]
        latents[i] = state.latents[i] # REVIEW: this vs `get_latents`
        kinematics[i] = get_state(elements[i], sim)
    end

    BulletState(elements, latents, kinematics)

end

end # module PhyBullet
