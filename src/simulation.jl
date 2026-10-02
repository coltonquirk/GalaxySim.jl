using LinearAlgebra

@doc raw"""
Particle type to hold position, velocity and mass information
"""
Base.@kwdef mutable struct Particle{
    Tv <: AbstractVector{<:AbstractFloat},
    Tm <: AbstractFloat
}
    pos :: Tv
    vel :: Tv
    m :: Tm
end

"""
Particles type to hold collection of all particles in the simulation
"""
Base.@kwdef mutable struct Particles{
  Tp <: AbstractVector{Particle}
}
    particles::Tp
end

@doc raw"""
grav_force(p1::Particle, p2::Particle, G::Float64=1.0)

Function to calculate the force of particle 1 from particle 2
F = -G M₁ M₂ / r^2
By default uses normalized units with G=1.0


# Examples
```jldoctest
julia> grav_force(p1, p2)
insert result of grav_force(particles)
```
"""
function grav_force(p1::Particle, p2::Particle, G::Float64=1.0)
    # Vector pointing from particle 2 to particle 1
    r = p1.pos - p2.pos
    force = -G * p1.m * p2.m * r / norm(r)^3
    return force
end

@doc raw"""
sum_forces(particles::Particles)

Function that sums the forces acting on the particles,
this can be from just the other particles interactinf gravitationally
or with additional force components such as a galactic potential.

# Examples
```jldoctest
julia> sum_forces(particles)
insert result of sum_forces(particles)
```
"""
function sum_forces(particles::Particles)
    grav_forces = fum_grav_forces(particles::Particles)
end

@doc raw"""
sum_grav_forces(particles::Particles)

Function that sums the gravitational forces from
all the particles
This can be done with:
- Direct summation O(N^2)
- Barnes-hut approximation ~O(N log(N))
- Fast Multipole Method O(N)

# Examples
```jldoctest
julia> sum_grav_forces(particles)
insert result of sum_grav_forces(particles)
```
"""
function sum_grav_forces(particles::Particles)

end

@doc raw"""
integrate(
    particles::Particles,
    forces::Vector{Vector{Float64}},
    dt::Float64
)

Function that takes the particles and integrates the forces applied
to the particles to update their position and velocities.
various integrators are implemented:
- euler
- 2nd order euler
- velocity-verlet/leapfrog
- rk2
- rk4
- yoshida

# Examples
```jldoctest
julia> integrate(particles, forces, dt)
insert result of integrate(particles, forces, dt)
```
"""
function integrate(
    particles::Particles,
    forces::Vector{Vector{Float64}},
    dt::Float64
)

end


@doc raw"""
update!(particles::Particles)

Updates the particles in the simulation in place.
Sums the forces acting on the particles.
Then implements an integration scheme to update the
position and velocity of the particles.

# Examples
```jldoctest
julia> dfun()
insert result of dfun()
```
"""
function update!(Particles)
    forces = sum_forces(particles)
    integrate(particles, forces)
end

@doc raw"""
seedparticles(N)

seed and N-body sim with N particles.

# Examples
```jldoctest
julia> seedparticles(100)
insert result of seedparticles(100)
```
"""
function seedparticles(N)

end

@doc raw"""
simulate(N::Int, t_start, t_end, dt)

Function that simulates an N-body problem with N particles
starting from t_start, ending at t_end, with timestep dt

# Examples
```jldoctest
julia> simulate()
insert result of simulate()
```
"""
function simulate(N::Int, t_start, t_end, dt)
    particles = seedparticles(N)

    t = t_start
    while t < t_end
        update!(particles, dt)
        t += dt
    end
end
