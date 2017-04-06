import rebound
import reboundx
import numpy as np
from collections import OrderedDict
from numpy.random import seed, normal, uniform

import sys
sys.path.append('../')

from ic import wrap, output, plot, removedamping, integrate, initialize, drawnormal, res_chain_setup

simID=int(sys.argv[1]) 
seed(simID)
K = 10**uniform(1,3)
  
incdist = {'b':[89.65,0.27], 'c':[89.67, 0.17], 'd':[89.75,0.16], 'e':[89.86, 0.12], 'f':[89.68,0.034], 'g':[89.710,0.025], 'h':[89.8, 0.3]}
for inc in incdist.values():
    inc[0] = (90-inc[0])*np.pi/180.
    inc[1] *= np.pi/180.

    massdist = {'b':[0.85,0.72], 'c':[1.38,0.61], 'd':[0.41,0.27], 'e':[0.62,0.58], 'f':[0.68,0.18], 'g':[1.34,0.88], 'h':[0.41, 0.27]}
Mearth = 3.e-6
Mstar = 0.08
for mass in massdist.values():
    mass[0] *= Mearth/Mstar
    mass[1] *= Mearth/Mstar

planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
resonances = OrderedDict([(('b','c'),(8,5)),(('c','d'),(5,3)),(('d','e'),(3,2)),(('e','f'),(3,2)),(('f','g'),(4,3)),(('g','h'),(3,2))]) # ordered so we add planets in right sequence
threebodyresonances = OrderedDict([(('b','c','d'),(2,3)),(('c','d','e'),(1,2)),(('d','e','f'),(2,3)),(('e','f','g'),(1,2)),(('f','g','h'),(1,1))])

sim = rebound.Simulation()
ps = sim.particles

rebx = reboundx.Extras(sim)
params = rebx.add("modify_orbits_forces")

delta=2.e-2 # fractional distance to start beyond resonance. Would be better to calculate relative to width of resonance

res_chain_setup(sim, planets, resonances, delta, massdist, incdist)

sim.integrator="whfast"
sim.dt=ps[1].P*0.05

outputs = initialize(planets, resonances, threebodyresonances)

taua = 5.e7 # 5e7/2pi = 8e6 inner planet orbits
taue = taua/K

for p in ps[1:]:
    p.params["tau_e"] = -taue
ps['h'].params["tau_a"] = -taua

T = taua*delta*50
Nout = 1000

filename = 'data/IC{0}K{1:.4e}'.format(simID, K)
sim.initSimulationArchive(filename+'.bin', interval=T/200)
times = np.linspace(0,T,Nout)
for i, time in enumerate(times):
    sim.integrate(time)
