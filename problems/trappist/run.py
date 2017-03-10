# python run.py 4res taue0 tmax where 4res means the starting binary is in data/4res.bin and taue0 is the initial taue for the innermost planet
import rebound
import reboundx
import numpy as np
import sys
import os

def collision(reb_sim, col):
    reb_sim.contents._status = 5
    return 0

simID = int(sys.argv[1])
tmax = float(sys.argv[2])

for filename in os.listdir('data/'):
    if filename.startswith("IC{0}K".format(simID)):
        path = 'data/'+filename
        sim = rebound.Simulation.from_archive(path)
        sim.simulationarchive_filename = path.encode('ascii')

sim.collision = "direct"
sim.collision_resolve = collision

sim.integrate(tmax)
