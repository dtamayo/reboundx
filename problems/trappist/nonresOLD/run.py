import rebound
import reboundx
import numpy as np
from collections import OrderedDict
from numpy.random import seed, normal, uniform

import sys
sys.path.append('../')
from ic import drawnormal

simID=int(sys.argv[1]) 
seed(simID)
  
incdist = {'b':[89.65,0.27], 'c':[89.67, 0.17], 'd':[89.75,0.16], 'e':[89.86, 0.12], 'f':[89.68,0.034], 'g':[89.710,0.025], 'h':[89.8, 0.3]}
for inc in incdist.values():
    inc[0] = (90-inc[0])*np.pi/180.
    inc[1] *= np.pi/180.

    massdist = {'b':[0.85,0.72], 'c':[1.38,0.61], 'd':[0.41,0.27], 'e':[0.62,0.58], 'f':[0.68,0.18], 'g':[1.34,0.88], 'h':[0.41, 0.27]}
Mearth = 3.e-6
Mstar = 0.08
for mass in massdist.values():
    mass[0] *= Mearth
    mass[1] *= Mearth

Pdist = {'b':[1.51087081,6.e-7], 'c':[2.4218233,1.7e-6], 'd':[4.049610,6.3e-5], 'e':[6.099615,1.1e-4], 'f':[9.206690,1.5e-4], 'g':[12.35294,1.2e-3], 'h':[18.764,9.e-3]} # days
for p in Pdist.values():
   p[0] /= 365.25
   p[1] /= 365.25

eccmax = {'b':0.081, 'c':0.083, 'd':0.07, 'e':0.085, 'f':0.063, 'g':0.061, 'h':0.07}

planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
resonances = OrderedDict([(('b','c'),(8,5)),(('c','d'),(5,3)),(('d','e'),(3,2)),(('e','f'),(3,2)),(('f','g'),(4,3)),(('g','h'),(3,2))]) # ordered so we add planets in right sequence
threebodyresonances = OrderedDict([(('b','c','d'),(2,3)),(('c','d','e'),(1,2)),(('d','e','f'),(2,3)),(('e','f','g'),(1,2)),(('f','g','h'),(1,1))])

sim = rebound.Simulation()
ps = sim.particles

sim.add(m=Mstar)
masses, incs, Omegas, thetas, periods, pomegas, eccs = {}, {}, {}, {}, {}, {}, {}

for p in planets:
    masses[p] = drawnormal(massdist[p])
    incs[p] = drawnormal(incdist[p])
    Omegas[p] = uniform(0,2*np.pi)
    thetas[p] = uniform(0,2*np.pi)
    eccs[p] = uniform(0,eccmax[p])
    pomegas[p] = uniform(0,2*np.pi)
    periods[p] = drawnormal(Pdist[p])
    sim.add(m=masses[p], P=periods[p], e=eccs[p], inc=incs[p], Omega=Omegas[p], pomega=pomegas[p], theta=thetas[p],hash=p)

sim.move_to_com()
ps = sim.particles
Hillradius = (ps[1].m/(3.*ps[0].m))**(1./3.)*ps[1].a
for p in sim.particles[1:]:
    p.r = Hillradius

def collision(reb_sim, col):
    reb_sim.contents._status = 5
    return 0
sim.collision = "direct"
sim.collision_resolve = collision
maxdistance=1. # AU
sim.exit_max_distance=maxdistance

sim.integrator="whfast"
sim.dt=ps[1].P*0.07

filename = 'data/nonresIC{0}'.format(simID)
sim.initSimulationArchive(filename+'.bin', interval=1000.)
tmax=5.e6
try:
    sim.integrate(tmax)
except:
    sim._status = 5
