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
sim.dt = ps[1].P*0.07/5.

aref = ps[1].a
mref = ps[1].m
for p in ps[1:]:
    p.params["tides_synchronous_tau_e"] = -taue*(p.m/mref)*(p.a/aref)**6.5

rebx.integrator="rk4"
damping.force_as_operator = 1
damping.operator_order = 1

filename='data/mag{0:.4e}taue{1:.2e}'.format(mag,taue)
rebx.save(filename+"rebx.bin")

tmax = 1.e6*(taue/100)/2.
interval = tmax/1000
sim.initSimulationArchive(filename+'.bin', interval=interval)
sim.integrate(tmax)
