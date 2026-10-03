using GalaxySim
using GLMakie

x = 0:1e-2:2π
y = sin.(x)

fig = Figure()
ax = Axis(fig[1,1])
lines!(x, y)

test1 = Particle(pos=[0.0, 0.0], vel=[0.0,0.0], m=1.0)
test2 = Particle(pos=[1.0, 0.0], vel=[0.0,1.5], m=1.0)

particles = Particles(particles = [test1, test2])

t_end = 2.0
euler_test(test2, test1; t_end=t_end)
