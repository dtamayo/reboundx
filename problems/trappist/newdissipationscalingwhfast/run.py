import rebound
import reboundx
import numpy as np
from collections import OrderedDict
from numpy.random import seed, uniform, normal
import sys

def wrap(val):
    while val < 0:
        val += 2*np.pi
    while val > 2*np.pi:
        val -= 2*np.pi
    return val*180/np.pi

taue = float(sys.argv[1])
mag = float(sys.argv[2])

filename='data/startermag{0:.4e}'.format(mag)
sim = rebound.Simulation.from_file(filename)

rebx = reboundx.Extras(sim)
damping = rebx.add("tides_synchronous_ecc_damping")

ps = sim.particles
aref = ps[1].a
mref = ps[1].m
for p in ps[1:]:
    p.params["tides_synchronous_tau_e"] = -taue*(p.m/mref)*(p.a/aref)**6.5

filename='data/mag{0:.4e}taue{1:.2e}'.format(mag,taue)
rebx.save(filename+"rebx.bin")

tmax = 5.e5*(taue/100)
interval = tmax/1000
sim.initSimulationArchive(filename+'.bin', interval=interval)
sim.integrate(tmax)
