# Generates a scene containing an elastic cube within a larger cubic enclosure
# Given its elasticity, geometry, and mass, etc. how will the cube bounce?
# Implemented as a deterministic simulation over uncertain properties
# Two ways to provide (fix) observed values - parameters and constraints

using Accessors
using Distributions
using Plots
using Gen
using PyCall
using PhySMC
using PhyBullet
bullet = pyimport("pybullet")
pybullet_data = pyimport("pybullet_data")
include("truncatednorm.jl")

# TODO: Test cases - make sure it can recover elasticity
# TODO: Find some way to sample from particles and compare


# Plots

@userplot SimPlot
@recipe function f(cp::SimPlot)
    z, t = cp.args
    cs = size(z, 1)
    k = 10
    inds = (max(1, t-k):t)
    n = length(inds)
    linewidth -->range(0,10, length = n)
    seriesalpha --> range(0,1,length = n)
    xguide --> "time"
    yguide --> "height of cube (z)"
    ylims --> (0, 4.0)
    xlims --> (1, 60)
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

function animate_trace(traces::Gen.Trace; label = "trace")
    t = first(get_args(trace))
    zs = reshape(get_zs(trace), (t, 1))
    anim = @animate for i = 2:t
        simplot(zs, i, label = label)
    end
end

function animate_traces(traces::Vector{<:Gen.Trace})
    n = length(traces)
    zzs = reduce(hcat, map(get_zs, traces))
    t = size(zzs, 1)
    anim = @animate for i=2:t
        simplot(zzs, i)
    end
end


# Generative Model

# Generate initial scene configuration
function init_scene(mass::Float64=1.0, restitution::Float64=0.9)

    bullet.setGravity(0, 0, -10)

    # Create and position cube
    cubeBody = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.2, .2, .2])
    startOrientationCube = bullet.getQuaternionFromEuler([0, 0, 1])
    cube = bullet.createMultiBody(baseCollisionShapeIndex=cubeBody, basePosition=[0., 0., 3.], baseOrientation=startOrientationCube)
    bullet.changeDynamics(cube, -1, mass=mass, restitution=restitution)

    # Create and position ground plane
    planeID = bullet.createCollisionShape(bullet.GEOM_PLANE)
    plane = bullet.createMultiBody(baseCollisionShapeIndex=planeID, basePosition=[0, 0, 0])
    bullet.changeDynamics(plane, -1, mass=0.0, restitution=0.9)

    # Create and position walls
    wall1ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.02, 1, 1])
    quaternion = bullet.getQuaternionFromEuler([0, 0, 0])
    plane2 = bullet.createMultiBody(baseCollisionShapeIndex=wall1ID, basePosition=[-1, 0, 1], baseOrientation=quaternion)
    bullet.changeDynamics(plane2, -1, mass=0.0, restitution=0.9)

    wall2ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[1, .02, 1])
    plane3 = bullet.createMultiBody(baseCollisionShapeIndex=wall2ID, basePosition=[0, 1, 1], baseOrientation=quaternion)
    bullet.changeDynamics(plane3, -1, mass=0.0, restitution=0.9)

    wall3ID = bullet.createCollisionShape(bullet.GEOM_BOX, halfExtents=[.02, 1, 1])
    plane4 = bullet.createMultiBody(baseCollisionShapeIndex=wall3ID, basePosition=[1, 0, 1], baseOrientation=quaternion)
    bullet.changeDynamics(plane4, -1, mass=0.0, restitution=0.9)

    return cube
end


# Update latent variables 
function update_latents(latents::RigidBodyLatents, mass::Float64, res::Float64)
    RigidBodyLatents(setproperties(latents.data, mass=mass, restitution=res))
end


# Sample initial estimates of mass and restitution from their priors
@gen function sample_from_prior(latents::RigidBodyLatents)
    mass = {:mass} ~ gamma(1.2, 10.)
    res = {:restitution} ~ uniform(0, 1)
    new_latents = update_latents(latents, mass, res)
    return new_latents
end


# Extract rigid body position and add noise
@gen function generate_observation(k::RigidBodyState)
    pos = k.position
    obs = {:position} ~ broadcasted_normal(pos, 0.1)
    return obs
end

# Deterministically generate the next state, then sample an observation given that state
@gen function kernel(t::Int, prev_state::BulletState, sim::BulletSim)
    next_state::BulletState = PhySMC.step(sim, prev_state)
    {:observation} ~ Gen.Map(generate_observation)(next_state.kinematics)
    return next_state
end

# Sample latents and run complete forward simulation
@gen function simulation(T::Int, sim::BulletSim, init_state::BulletState)
    latents = {:prior} ~ Gen.Map(sample_from_prior)(init_state.latents)
    init_state = Accessors.setproperties(init_state; latents=latents)
    states = {:kernel} ~ Gen.Unfold(kernel)(T, init_state, sim)
    return states
end


# Inference Procedure

# Rejuvenate latent estimates
@gen function proposal(trace::Gen.Trace)

    # Parse trace for previous mass and restitution estimates
    choices  = get_choices(trace)
    prev_mass = choices[:prior => 1 => :mass]
    prev_res  = choices[:prior => 1 => :restitution]

    # Resample mass and restitution from Gaussians centered at their previous values
    mass = {:prior => 1 => :mass} ~ trunc_norm(prev_mass, 0.1, 0.0, Inf)
    restitution = {:prior => 1 => :restitution} ~ trunc_norm(prev_res, 0.1, 0.0, 1.0)

    return (mass, restitution)
end

# Particle filter
function infer(gm_args::Tuple, obs::Vector{Gen.ChoiceMap}, num_particles::Int=20)

    get_args(t) = (t, gm_args[2:3]...)
    state = Gen.initialize_particle_filter(simulation, get_args(0), EmptyChoiceMap(), num_particles)
    argdiffs = (UnknownChange(), NoChange(), NoChange())

    for (t, obs) = enumerate(obs)
        @elapsed begin
            for i=1:num_particles
                state.traces[i], _ = Gen.mh(state.traces[i], proposal, ())
            end
            
            Gen.maybe_resample!(state, ess_threshold=num_particles/2)
            Gen.particle_filter_step!(state, get_args(t), argdiffs, obs)
        end
    end

    return state.traces
end


# Main Function

function main()
    client = bullet.connect(bullet.DIRECT)::Int64
    bullet.setAdditionalSearchPath(pybullet_data.getDataPath())
    bullet.resetDebugVisualizerCamera(5, 0, -5, [0, 0, 2])

    # Initialize scene and simulation context
    rb_cube = RigidBody(init_scene())
    sim = BulletSim(; client=client)
    init_state = BulletState(sim, [rb_cube])

    # Generate ground truth trajectory
    args = (60, sim, init_state)
    gt_constraints = choicemap((:prior => 1 => :restitution, 0.8), (:prior => 1 => :mass, 1.0))
    ground_truth = first(generate(simulation, args, gt_constraints))

    # Read position observations from trace and store in a vector
    gt_choices = get_choices(ground_truth)
    t = args[1]
    observations = Vector{Gen.ChoiceMap}(undef, t)
    for i = 1:t
        obs = choicemap()
        address = :kernel => i => :observe
        set_submap!(obs, address, get_submap(gt_choices, address))
        observations[i] = obs
    end

    # Execute inference
    result = infer(args, observations)

    # Visualize particles
    gif(animate_traces(result), fps=24)


    # bullet.resetSimulation(bullet.RESET_USE_DEFORMABLE_WORLD)

    # gif(animate_trace(gt), fps=24)

    # run_sample(sphereBodies, sphere_locations)

    # traces = [first(generate(simulate, gargs)) for _=1:5]
    # anim = animate_traces(traces)
    # gif(anim, fps = 24)

    bullet.disconnect()
end


main()

#=

cube_location, platform_location, platform_orientation = sample_init_state()
cubeBodies = []

cube_locations = Vector{Float64}[]
push!(cube_locations, [0., 0., 3.])

push!(cubeBodies, cubeBody)
startPosCube = cube_locations[1]

cube1 = p.loadSoftBody("soft_cube__sf.obj", simFileName="soft_cube.vtk", basePosition=cube_locations[0],
                        mass=4, useNeoHookean=1, NeoHookeanMu=60, NeoHookeanLambda=600, NeoHookeanDamping=0.001,
                        collisionMargin=0.006, useSelfCollision=1, frictionCoeff=0.5, repulsionStiffness=800)
cubeBodies.append(cube1)

cube2 = p.loadSoftBody("soft_cube__sf.obj", simFileName="soft_cube.vtk", basePosition=cube_locations[1],
                        mass=4, useNeoHookean=1, NeoHookeanMu=280, NeoHookeanLambda=600, NeoHookeanDamping=0.001,
                        collisionMargin=0.006, useSelfCollision=1, frictionCoeff=0.5, repulsionStiffness=700)
cubeBodies.append(cube2)

# Create and position cube
cubeBody1 = p.createCollisionShape(p.GEOM_BOX, halfExtents=[.1, .1, .1])
cubeBodies.append(cubeBody1)
startPosCube = cube_locations[1]
startOrientationCube = p.getQuaternionFromEuler([0, 0, 1])
cube1 = p.createMultiBody(baseCollisionShapeIndex=cubeBody1, basePosition=startPosCube,
                          baseOrientation=startOrientationCube)
p.changeDynamics(cube1, -1, mass=1.0, restitution=0.7)

# Create and position cube
cubeBody2 = p.createCollisionShape(p.GEOM_BOX, halfExtents=[.1, .1, .1])
cubeBodies.append(cubeBody2)
startPosCube = cube_locations[2]
startOrientationCube = p.getQuaternionFromEuler([0, 0, 1])
cube2 = p.createMultiBody(baseCollisionShapeIndex=cubeBody2, basePosition=startPosCube,
                          baseOrientation=startOrientationCube)
p.changeDynamics(cube2, -1, mass=1.0, restitution=0.4)

# Create and position cube
cubeBody3 = p.createCollisionShape(p.GEOM_BOX, halfExtents=[.1, .1, .1])
cubeBodies.append(cubeBody3)
startPosCube = cube_locations[3]
startOrientationCube = p.getQuaternionFromEuler([0, 0, 1])
cube3 = p.createMultiBody(baseCollisionShapeIndex=cubeBody3, basePosition=startPosCube,
                          baseOrientation=startOrientationCube)
p.changeDynamics(cube3, -1, mass=1.0, restitution=0.1)

# Create and position platform
platformBody = p.createCollisionShape(p.GEOM_BOX, halfExtents=[.1, .1, .1])
startPosPlatform = platform_location
startOrientationPlatform = p.getQuaternionFromEuler([platform_orientation[0], platform_orientation[1],
                                                     platform_orientation[2]])
platform = p.createMultiBody(baseCollisionShapeIndex=platformBody, basePosition=startPosPlatform,
                             baseOrientation=startOrientationPlatform)
p.changeDynamics(platform, -1, mass=0.0, restitution=0.9)
=#

#=
function sample_init_state()
    # Generate platform's location from a range / uniform distribution
    random.seed()
    x = random.uniform(-2, 2)
    z = random.uniform(-2, 2)
    platform_location = [x, z, 2]

    # Condition cube's location on platform location
    x = random.uniform(x - .8, x + .8)
    z = random.uniform(z - .3, z + .3)
    cube_location = [x, z, 6]

    # Generate platform orientation (Von Mises)
    theta_y = random.vonmisesvariate(0, 30.06)
    platform_orientation = [0, theta_y, 0]

    return cube_location, platform_location, platform_orientation

end
=#

#=
function run_sample(cubeBodies, cube_locations)

    cube_positions = cube_locations

    for i=1:1000
        bullet.stepSimulation()
        current_position = bullet.getBasePositionAndOrientation(cubeBodies[1])[1]
        position_array = [i for i in current_position]
        push!(cube_positions, position_array)
        sleep(.001)
    end

    #trace = {"cube_position": cube_positions}

    #return trace
end
=#
