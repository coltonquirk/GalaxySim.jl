using GLMakie
using LinearAlgebra
using CSV
using Unitful, UnitfulAstro

const G = 1.0

mutable struct Particle
    pos
    vel
    m
end

function grav_force(p::Particle, p_other::Particle)
    r = p.pos - p_other.pos
    F = - G * p.m * p_other.m * r / norm(r)^3
    return F
end

function calc_forces(particles, force)
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

function integrate!(t, particles, force, dt)
    fR = calc_forces(particles, force)

    # Kick
    for (i, p) in enumerate(particles)
        p.vel += 0.5 * fR[i, :] * dt
    end

    # Drift
    for (i, p) in enumerate(particles)
        p.pos += p.vel * dt
    end

    fR_new = calc_forces(particles, force)

    # Kick
    for (i, p) in enumerate(particles)
        p.vel += 0.5 * fR_new[i, :] * dt
    end
end

function simulate(particles, force; t0=0.0, tf=20*4π, Δt=1e-2)
    t = t0
    hist = []
    while t<tf
        integrate!(t, particles, force, Δt)
        push!(hist, [particles[1].pos; particles[2].pos])
        t += Δt
    end

    reduce(hcat, hist)'
end

function plot_hist(hist)
    fig = Figure()
    ax = Axis3(fig[1,1])

    traj_1 = hist[:, 1:3] #  Point3f.(hist[:, 1:3])
    traj_2 = hist[:, 4:6] # Point3f.(hist[:, 4:6])

    lines!(ax, traj_1)
    lines!(ax, traj_2)

    fig
end

p1 = Particle(
    [0.0, 0, 0],
    [0, 0.0, 0],
    1.0
)

p2 = Particle(
    [1, 0, 0],
    [0, 2π, 0],
    1.0e-8
)

particles = [p1, p2]

test = simulate(particles, grav_force;)
fig = plot_hist(test)
