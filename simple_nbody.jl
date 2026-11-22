using GLMakie
using LinearAlgebra

const G = 1.0

mutable struct Particle
    pos
    vel
    m
end

mutable struct Particles{
    tv<:AbstractVector{Real},
    tm<:AbstractVector{Real}
}
    pos::tv
    vel::tv
    m::tm
end

mutable struct Particles{ps::Vector{Particle}}
    pos
    vel
    m
    function Particles(ps::Vecor{Particle})
        positions = zeros(length(ps), 3)
        velocities = zeros(length(ps), 3)
        masses = zeros(length(ps))

        for (i, p) in enumerate(ps)
            positions[i, :] = p.pos
            velocities[i, :] = p.vel
            masses[i] = p.m
        end

        new(positions, velocities, masses)
    end
end

function grav_force(p, p_other)
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
                forces[i, :] += f
                forces[j, :] -= f
            else
                continue
            end
        end
    end

    return forces
end

# TODO: Particles type
# particles getindex method
# Make sure to pass through the masses of the particles
# Figure out how to update the particles in place
function integrate!(t, particles, force, dt)
    fR = calc_forces(particles, force)

    # Kick
    # vel_half = zeros(length(particles), 3)
    for (i, p) in enumerate(particles)
        p.vel = @. p.vel + 0.5 * fR[i] * dt
    end

    # Drift
    # pos_new = zeros(length(particles), 3)
    for (i, p) in enumerate(particles)
        p.pos = @. p.pos + p.vel * dt
    end

    # Kick
    # fR_new = force(t, Particles(pos_new, vel_half, particles.m))
    fR_new = calc_forces(particles, force)
    # vel_new = zeros(length(particles), 3)
    for (i, p) in enumerate(particles)
        # vel_new[i, :] = vel_half[i, :] + 0.5 * fR_new[i, :] * dt
        p.vel = @. p.vel + 0.5 * fR_new[i, :] * dt
    end

    # return Particles(pos_new, vel_new, particles.m)
end

p1 = Particle(
    [-1, 0, 0],
    [0, -1, 0],
    1.0
)

p2 = Particle(
    [1, 0, 0],
    [0, 1, 0],
    1.0
)

p3 = Particle(
    [0, 1, 0],
    [-1, 0, 0],
    1.0
)

# force(p1, p2)
# force(p2, p1)

particles = [p1, p2, p3]

# [[0, 0, 0] for n in 1:length(particles)]
# test = zeros((length(particles), 3))
# test[1, :] += [1, 2, 3]
# test
# calc_forces(particles)

integrate!(0.0, particles, grav_force, 1e-2)
particles

t0 = 0.0
tf = 20.0
t = t0
Δt = 1e-2

while t<tf


    t += Δt
end
