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

###############################################################################
#                                   Plotting                                  #
###############################################################################


function plot_hist(hist; skip=100)
    N = size(hist)[2] # ÷3
    fig = Figure()
    ax = Axis3(fig[1,1])

    for i in 1:3:(N-2)
        lines!(ax, hist[begin:skip:end, i:i+2])
    end

    fig
end

function anim_hist(hist; limits=(nothing, nothing, nothing), skip=100, seconds=1e-2)
    idx = Observable(1)

    fig = Figure()
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
