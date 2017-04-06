import rebound
import reboundx
import numpy as np
import sys
import os
import re

def collision(reb_sim, col):
    reb_sim.contents._status = 5
    return 0

simID = int(sys.argv[1])
tmax = float(sys.argv[2])

path = 'data/7pdissipation/'

for filename in os.listdir(path):
    result = re.search("IC{0}K(.*)[0-9].bin".format(simID), filename)
    if result:
        sim = rebound.Simulation.from_archive(path+filename)
        rebx = reboundx.Extras.from_file(sim, path+filename[:-4]+'rebx.bin')
        damping = rebx.get_effect("tides_synchronous_ecc_damping")
        rebx.integrator="rk4"
        damping.force_as_operator = 1
        break

sim.collision = "direct"
sim.collision_resolve = collision

sim.integrate(tmax)
