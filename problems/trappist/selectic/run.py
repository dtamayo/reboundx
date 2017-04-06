import rebound
import reboundx
import numpy as np
import os
import sys
import pandas as pd
import re
from numpy.random import seed, normal, uniform
from numpy.testing import assert_approx_equal

sys.path.append('../')

from ic import output, initialize, removedamping, integrate

from collections import OrderedDict
planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
resonances = OrderedDict([(('b','c'),(8,5)),(('c','d'),(5,3)),(('d','e'),(3,2)),(('e','f'),(3,2)),(('f','g'),(4,3)),(('g','h'),(3,2))]) # ordered so we add planets in right sequence
threebodyresonances = OrderedDict([(('b','c','d'),(2,3)),(('c','d','e'),(1,2)),(('d','e','f'),(2,3)),(('e','f','g'),(1,2)),(('f','g','h'),(1,1))])

dfID = int(sys.argv[1])
idstart = int(sys.argv[2])
idend = int(sys.argv[3])
simIDs = range(idstart, idend)

df = pd.DataFrame(columns=('K','mag','filename'))

frac = 0.2 # take last x fraction of samples to check for resonance
pratiothresh = 0.02

path = '../makeic/data/'

for simID in simIDs:
    outputs = initialize(planets, resonances, threebodyresonances)
    for filename in os.listdir(path):
        result = re.search('IC(.*)K(.*).bin', filename)
        if result:
            ID = int(result.group(1))
            if ID == simID:
                K = float(result.group(2))
                sa = rebound.SimulationArchive(path+filename)
                Nout = len(sa)
                N = sa[0].N
                for i,sim in enumerate(sa):
                    ps = sim.particles
                    output(sim,planets,resonances,threebodyresonances,outputs)

                t, e, P, Pratio, phi1, phi2, deltapomega, phi3body = outputs

                success = True
                for resonance in resonances.items():
                    pair = resonance[0]
                    res = resonance[1]
                    resratio = res[0]/res[1]
                    startindex = -int(len(Pratio[pair])*frac)
                    pratio = np.array(Pratio[pair][startindex:])
                    if abs((pratio.mean()-resratio)/resratio) > pratiothresh:
                        print("{0} pair not in {1} resonance. Mean period ratio = {2}".format(pair, res, pratio.mean()))
                        success = False

                if success:
                    sim = sa[-1]
                    ps = sim.particles

                    rebx = reboundx.Extras(sim)
                    params = rebx.add("modify_orbits_forces")

                    taua = 5.e7
                    taue = taua/K

                    for p in ps[1:]:
                        p.params["tau_e"] = -taue
                    ps['h'].params["tau_a"] = -taua

                    Tremoval = 5 # number of taue0s over which to remove damping
                    removedamping(sim, Tremoval, K, planets,resonances,threebodyresonances,outputs)

                    seed(simID)
                    assert_approx_equal(K, 10**uniform(1,3), 4, "simID did not generate same value of K as in ic") # get random number after K for seed
                    mag = 10**uniform(-3,0)

                    Mstar = 0.08
                    sim2 = rebound.Simulation()
                    sim2.G = 4*np.pi**2
                    sim2.add(m=Mstar)

                    Pd = 4.0496/365.25

                    mscale = Mstar/sim.particles[0].m
                    Pscale = Pd/sim.particles['d'].P

                    ps = sim.particles
                    for p in planets: 
                        r = uniform(-1,1)
                        sign = abs(r)/r # get +/- 1 randomly
                        sim2.add(m=ps[p].m*mscale, P=ps[p].P*Pscale, e=ps[p].e*(1. + sign*mag), inc=ps[p].inc, pomega=ps[p].pomega, Omega=ps[p].Omega, theta=ps[p].theta, hash=p)

                    sim2.move_to_com()
                    sim = sim2
                    ps = sim.particles
                    sim.integrator="whfast"
                    sim.dt=ps[1].P*0.07 

                    filename='IC{0}K{1:.4e}mag{2:.4e}.bin'.format(simID, K, mag)
                    df.loc[simID] = [K,mag,filename]

                    sim.initSimulationArchive('data/'+filename, interval=1.e3)
                    sim.integrate(2.e3, exact_finish_time=0)
                break

df.to_csv('data/df{0}.csv'.format(dfID))
