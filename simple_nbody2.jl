using GLMakie
using LinearAlgebra

const G = 1.0

mutable struct Particle{
    Tmat<:AbstractVecOrMat{<:Real},
    Tmass<:Real
}
    pos::Tmat
    vel::Tmat
    m::Tmass
end

function grav_force(p::Particle, p_other::Particle)
    r = p.pos - p_other.pos
    F = - G * p.m * p_other.m * r / norm(r)^3
    return F
end

function calc_forces(particles::Vector{<:Particle}, force::Function)
    forces = zeros((length(particles), 3))
    for (i, p) in enumerate(particles)
        for (j, p_other) in enumerate(particles)
            if i > j
                f = force(p, p_other)
                # Need to calculate acceleration not force
                forces[i, :] += f / p.m
                forces[j, :] -= f / p_other.m
            else
                continue
            end
        end
    end
    return forces
end

function integrate!(t::Float64, particles::Vector{<:Particle}, force::Function, Δt::Float64)
    fR = calc_forces(particles, force)

    # TODO: Can I combine the for loops? There is something I vaguely remember about different timesteps/step
    # Kick
    for (i, p) in enumerate(particles)
        p.vel += 0.5 * fR[i, :] * Δt
    end

    # Drift
    for (i, p) in enumerate(particles)
        p.pos += p.vel * Δt
    end

    fR_new = calc_forces(particles, force)
    # Kick
    for (i, p) in enumerate(particles)
        p.vel += 0.5 * fR_new[i, :] * Δt
    end
end

function simulate(
    particles::Vector{<:Particle},
    force::Function;
    t0::Float64=0.0,
    tf::Float64=20*4π,
    Δt::Float64=1e-2
)
    N = length(particles)
    t = t0
    hist::Vector{Vector{Float64}} = [Vector{Float64}(undef,N*3)]
    while t<tf
        integrate!(t, particles, force, Δt)
        push!(hist, vcat([p.pos for p in particles]...)) # alternative: push!(hist, [([p.pos for p in particles]...)...])
        t += Δt
    end
    reduce(hcat, hist[2:end])' # need to remove the random starting position # hist # stack(hist)
end

function simulate(
    orbit::Symbol=:default,
    force::Function=grav_force;
    t0::Float64=0.0,
    tf::Float64=20*4π,
    Δt::Float64=1e-2
)
    particles = initialize(orbit)
    simulate(particles, force; t0=t0, tf=tf, Δt=Δt)
end

###############################################################################
#                              Initial Conditions                             #
###############################################################################
using Random
rng = Random.default_rng()
Random.seed!(42)

function initialize(orbit::Symbol=:default)
    @assert in(orbit, [:default, :solar_system, :broucke, :figure_eight, :random]) "orbit: $(orbit) is not a valid condition"

    if orbit == :solar_system
        particles = solar_system()
    elseif orbit == :broucke
        particles = broucke()
    elseif orbit == :figure_eight
        particles = figure_eight()
    elseif orbit == :random
        particles = random_particles()
    else
        particles = default_particles()
    end

    return particles
end

function default_particles()
    p1 = Particle(
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        1.0
    )

    p2 = Particle(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        1.0
    )

    p3 = Particle(
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        1.0
    )

    return [p1, p2, p3]
end

using Unitful, UnitfulAstro
const M☼_kg = uconvert(u"kg", 1.0u"Msun")

function solar_system()
    Sun = Particle(
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        1.0
    )

    Mercury = Particle(
        [0.387, 0.0, 0.0],
        [0.0, √(1.0 / 0.387), 0.0], # Calculate circular velocity of an orbit. √(GM/r)
        ustrip(0.330e24u"kg"/M☼_kg)
    )

    Venus = Particle(
        [0.723, 0.0, 0.0],
        [0.0, √(1.0 / 0.723), 0.0], # Calculate circular velocity of an orbit. √(GM/r)
        ustrip(4.87e24u"kg"/M☼_kg)
    )

    Earth = Particle(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0], # Why does 1.0 work for the simulated Earth example? Shouldn't it be 2 or 4π?
        # Circular orbital velocity with G=M=1.0 (technically GM should be 4π²)
        # v = √(GM/r) = √(1*1/1)=1.0
        ustrip(5.97e24u"kg"/M☼_kg)
    )

    Mars = Particle(
        [1.524, 0.0, 0.0],
        [0.0, √(1.0 / 1.524), 0.0], # Calculate circular velocity of an orbit. √(GM/r)
        ustrip(0.642e24u"kg"/M☼_kg)
    )

    Jupiter = Particle(
        [5.2044, 0.0, 0.0],
        [0.0, √(1.0 / 5.2044), 0.0], # Calculate circular velocity of an orbit. √(GM/r)
        ustrip(1898.0e24u"kg"/M☼_kg)
    )

    Saturn = Particle(
        [9.5826, 0.0, 0.0],
        [0.0, √(1.0 / 9.5826), 0.0], # Calculate circular velocity of an orbit. √(GM/r)
        ustrip(1898.0e24u"kg"/M☼_kg)
    )

    Uranus = Particle(
        [19.21840, 0.0, 0.0],
        [0.0, √(1.0 / 19.21840), 0.0], # Calculate circular velocity of an orbit. √(GM/r)
        ustrip(86.8e24u"kg"/M☼_kg)
    )

    Neptune = Particle(
        [30.11, 0.0, 0.0],
        [0.0, √(1.0 / 30.11), 0.0], # Calculate circular velocity of an orbit. √(GM/r)
        ustrip(102.0e24u"kg"/M☼_kg)
    )

    return [Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune]
end

function figure_eight()
    # Orbit I.A.1 i.c.
    p1 = Particle(
        [-1.0, 0.0, 0.0],
        [0.3471168881, 0.5327249454, 0.0],
        1.0
    )

    p2 = Particle(
        [1.0, 0.0, 0.0],
        [0.3471168881, 0.5327249454, 0.0],
        1.0
    )

    p3 = Particle(
        [0.0, 0.0, 0.0],
        [-0.6942337762, -1.0654498908, 0.0],
        1.0
    )

    return [p1, p2, p3]
end

function broucke()
    # Broucke A 1
    p1 = Particle(
        [-0.9892620043, 0.0, 0.0],
        [0.0, 1.9169244185, 0.0],
        1.0
    )

    p2 = Particle(
        [2.2096177241, 0.0, 0.0],
        [0.0, 0.1910268738, 0.0],
        1.0
    )

    p3 = Particle(
        [-1.2203557197, 0.0, 0.0],
        [0.0, -2.1079512924, 0.0],
        1.0
    )

    return [p1, p2, p3]
end

function random_particles(N::Integer=3)
    particles::Vector{<:Particle} = [Particle([0.0,0.0,0.0], [0.0,0.0,0.0], 0.0)]
    for _ in 1:N
        p = Particle(
            2 .* rand(3) .- 1,
            2 .* rand(3) .- 1,
            rand()
        )
        push!(particles, p)
    end
    return particles[2:end]
end

###############################################################################
#                                   Plotting                                  #
###############################################################################
using DataStructures

const inch = 96
const cm = inch / 2.54
const pt = 4/3

black_latexfonts = merge(theme_black(), theme_latexfonts());
set_theme!(black_latexfonts)

function plot_hist(hist; skip=100)
    N = size(hist)[2] # ÷3
    fig = Figure(size=(7inch, 10inch), fontsize=24pt)
    ax = Axis3(fig[1,1])

    for i in 1:3:(N-2)
        lines!(ax, hist[begin:skip:end, i:i+2])
    end

    fig
end

function anim_hist(hist; limits=(nothing, nothing, nothing), skip=100, seconds=1e-2)
    idx = Observable(1)

    fig = Figure(size=(10inch, 7inch), fontsize=24pt)
    ax = Axis3(
        fig[1,1],
        limits = limits
    )

    data = @lift hist[$idx, :]

    points = @lift [Point3f($data[i:i+2]) for i in 1:3:(length($data)-2)]
    scatter!(points)

    display(fig)
    for i in 1:skip:size(hist)[1]
        idx[] = i
        # autolimits!(ax)
        sleep(seconds)
    end
end

###############################################################################
#                               Testing/Example                               #
###############################################################################

p1 = Particle(
    [-1.0, 0.0, 0.0],
    [0.0, -1.0, 0.0],
    1.0
)

p2 = Particle(
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    1.0
)

p3 = Particle(
    [0.0, 1.0, 0.0],
    [-1.0, 0.0, 0.0],
    1.0
)

particles = [p1, p2, p3]
test = simulate(particles, grav_force;)
plot_hist(test)
anim_hist(
    test;
    limits = (-90, 10, -10, 10, -2, 2),
    skip=50,
    seconds = 1e-2
)

Sun = Particle(
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    1.0
)

Earth = Particle(
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0], # Why does 1.0 work for the simulated Earth example? Shouldn't it be 2 or 4π?
    # TODO: Circular orbital velocity with G=M=1.0 (technically GM should be 4π)
    # v = √(GM/r) = √(1*1/1)=1.0
    1e-8
)

solar_system = [Sun, Earth]
sol = simulate(solar_system, grav_force;)
plot_hist(sol; skip=10)
anim_hist(
    sol;
    limits = (-1, 1, -1, 1, -1, 1),
    skip=10,
    seconds=1e-2
)


# Orbit I.A.1 i.c.
p1 = Particle(
    [-1.0, 0.0, 0.0],
    [0.3471168881, 0.5327249454, 0.0],
    1.0
)

p2 = Particle(
    [1.0, 0.0, 0.0],
    [0.3471168881, 0.5327249454, 0.0],
    1.0
)

p3 = Particle(
    [0.0, 0.0, 0.0],
    [-0.6942337762, -1.0654498908, 0.0],
    1.0
)

particles = [p1, p2, p3]
periodic = simulate(particles, grav_force)
plot_hist(periodic)

# Broucke A 1
p1 = Particle(
    [-0.9892620043, 0.0, 0.0],
    [0.0, 1.9169244185, 0.0],
    1.0
)

p2 = Particle(
    [2.2096177241, 0.0, 0.0],
    [0.0, 0.1910268738, 0.0],
    1.0
)

p3 = Particle(
    [-1.2203557197, 0.0, 0.0],
    [0.0, -2.1079512924, 0.0],
    1.0
)

particles = [p1, p2, p3]
broucke = simulate(particles, grav_force; tf = 20.0, Δt=1e-4)
plot_hist(broucke; skip=500)
anim_hist(
    broucke;
    limits = (-3, 3, -1, 1, -2, 2),
    skip=50,
    seconds = 1e-4
)
