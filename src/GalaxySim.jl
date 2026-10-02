module GalaxySim

using LinearAlgebra

include("simulation.jl")
include("barnes_hut/barnes_hut.jl")

export Particle
export Particles
export grav_force
export euler_test

end # module GalaxySim
