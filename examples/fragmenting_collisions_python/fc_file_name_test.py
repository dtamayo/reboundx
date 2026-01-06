import rebound
import reboundx
import unittest
from ctypes import c_char_p

# Merging
sim = rebound.Simulation()
rebx = reboundx.Extras(sim)
collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
rebx.add_collision_resolve(collision_resolve)
collision_resolve.params["fc_min_frag_mass"] = 1.0
# Write python string to param
collision_resolve.params["fc_particle_list_file"] = "family_tree.csv"
# Read param and convert to python string
print("Read: ", collision_resolve.params["fc_particle_list_file"])

sim.add(m=1, r=1)
sim.add(m=1, r=1,x=2.5,vx=-1)
sim.dt = 1
sim.integrator = "mercurius"
sim.collision = "direct"

print("Before collision, N = ", sim.N)

sim.step()

print("After collision, N = ", sim.N)

