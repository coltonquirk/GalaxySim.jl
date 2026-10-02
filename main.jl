using GalaxySim
using GLMakie

test1 = Particle(pos=[0.0, 0.0], vel=[0.0,0.0], m=1.0)
test2 = Particle(pos=[1.0, 0.0], vel=[0.0,0.0], m=1.0)
grav_force(test1, test2)

dt = 1e-2
test1.vel += grav_force(test1, test2)*dt
test1.pos += test1.vel * dt
