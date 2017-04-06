import rebound
import reboundx
import numpy as np
import sys
import os
import re
import pandas as pd

def collision(reb_sim, col):
    reb_sim.contents._status = 5
    return 0

filename = sys.argv[1]
tmax = float(sys.argv[2])

path = '../selectic/data/' + filename
sim = rebound.Simulation.from_archive(path)
sim.simulationarchive_filename = path.encode('ascii')

sim.collision = "direct"
sim.collision_resolve = collision

sim.integrate(tmax)
