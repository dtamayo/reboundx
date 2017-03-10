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

def removedamping(sim, Tremoval, K): # Tremoval in # of taue0s
    ps = sim.particles
    taue0 = ps[-1].params["tau_e"]
    T0 = sim.t
    T = abs(taue0)*Tremoval
    Nout = 1000
    times = np.linspace(T0, T0+T, Nout)
    
    for i, time in enumerate(times):
        for p in ps[1:]:
            p.params["tau_e"] = taue0/(1.-(sim.t-T0)/T)
            try: # try except so we don't assign a tau_a to planets that don't have it
                has_tau_a = p.params["tau_a"]
                p.params["tau_a"] = taue0*K/(1.-(sim.t-T0)/T)
            except:
                pass
        sim.integrate(time)
        
    for p in ps[1:]:
        p.params["tau_e"] = np.inf
        p.params["tau_a"] = np.inf

def drawnormal(pair):
    mean, sigma = pair
    val = -1
    while val < 0:
        val = normal(mean, sigma)
    return val

def res_chain_setup(sim, planets, resonances, delta):
    masses, incs, Omegas, thetas = {}, {}, {}, {}

    for p in planets:
        masses[p] = drawnormal(massdist[p])
        incs[p] = drawnormal(incdist[p])
        Omegas[p] = uniform(0,2*np.pi)
        thetas[p] = uniform(0,2*np.pi)

    ps = sim.particles
    sim.add(m=1.)
    
    p0 = planets[0]
    sim.add(m=masses[p0],a=1., inc=incs[p0], Omega=Omegas[p0], theta=thetas[p0], hash=p0)
    
    for resonance in resonances.items():
        pair = resonance[0]
        p1 = pair[0]
        p2 = pair[1]
        res = resonance[1]
        p = res[1]
        q = res[0]-res[1]
        sim.add(m=masses[p2],P=(p+q)/p*sim.particles[-1].P*(1.+delta), inc=incs[p2], Omega=Omegas[p2], theta=thetas[p2], hash=p2)
        
    sim.move_to_com() # Moves to the center of momentum frame

incdist = {'b':[89.65,0.27], 'c':[89.67, 0.17], 'd':[89.75,0.16], 'e':[89.86, 0.12], 'f':[89.68,0.034], 'g':[89.710,0.025]}
for inc in incdist.values():
    inc[0] = (90-inc[0])*np.pi/180.
    inc[1] *= np.pi/180.

    massdist = {'b':[0.85,0.72], 'c':[1.38,0.61], 'd':[0.41,0.27], 'e':[0.62,0.58], 'f':[0.68,0.18], 'g':[1.34,0.88]}
Mearth = 3.e-6
Mstar = 0.08
for mass in massdist.values():
    mass[0] *= Mearth/Mstar
    mass[1] *= Mearth/Mstar

planets = ['d', 'e', 'f', 'g']
resonances = OrderedDict([(('d','e'),(3,2)),(('e','f'),(3,2)),(('f','g'),(4,3))]) # ordered so we add planets in right sequence

simID=int(sys.argv[1])
seed(simID)
K = 10**uniform(1,3)
sim = rebound.Simulation()
ps = sim.particles

rebx = reboundx.Extras(sim)
params = rebx.add("modify_orbits_forces")

delta=2.e-2 # fractional distance to start beyond resonance. Would be better to calculate relative to width of resonance

res_chain_setup(sim, planets, resonances, delta)

sim.integrator="whfast"
sim.dt=ps[1].P/20.

taua = 1.e7
taue = taua/K

for p in ps[1:]:
    p.params["tau_e"] = -taue
ps['g'].params["tau_a"] = -taua

T = taua*delta*30
Nout = 1000

times = np.linspace(sim.t,sim.t+T,Nout)
for i, time in enumerate(times):
    sim.integrate(time)

Tremoval = 5 # number of taue0s over which to remove damping
removedamping(sim, Tremoval, K)

Pb = 1.511/365.25
Pc = 2.4218/365.25
Pg = 12.353/365.25

mscale = Mstar/sim.particles[0].m
Pscale = Pg/sim.particles[-1].P

sim2 = rebound.Simulation()
sim2.G = 4*np.pi**2
sim2.add(m=Mstar)
sim2.add(m=1.e-9, P=Pb, inc=drawnormal(incdist['b']), Omega = uniform(0,2*np.pi), theta=uniform(0,2*np.pi), hash='b')
sim2.add(m=1.e-9, P=Pc, inc=drawnormal(incdist['c']), Omega = uniform(0,2*np.pi), theta=uniform(0,2*np.pi), hash='c')

for p in planets:
    p1 = sim.particles[p]
    sim2.add(m=p1.m*mscale, P=p1.P*Pscale, e=p1.e, inc=p1.inc, pomega=p1.pomega, Omega=p1.Omega, theta=p1.theta, hash=p)

sim = sim2
ps = sim.particles
sim.integrator="whfast"
sim.dt=ps[1].P*0.07

planets = ['b', 'c', 'd', 'e', 'f', 'g']
resonances = OrderedDict([(('b','c'),(8,5)), (('c','d'),(5,3)), (('d','e'),(3,2)),(('e','f'),(3,2)),(('f','g'),(4,3))]) # ordered so we add planets in right sequence

T0 = sim.t
T = 1e5*ps[1].P
Nout = 1000
times = np.linspace(T0, T0+T, Nout)
        
Mb = drawnormal(massdist['b'])*mscale
Mc = drawnormal(massdist['c'])*mscale
 
for i, time in enumerate(times):
    ps['b'].m = Mb*(sim.t-T0)/T
    ps['c'].m = Mc*(sim.t-T0)/T
    sim.integrate(time)

mag =10**uniform(-3,0)
for p in planets: 
    theta = uniform(0,2*np.pi)
    ps[p].vx += mag*np.cos(theta)*ps[p].e*ps[p].v
    ps[p].vy += mag*np.sin(theta)*ps[p].e*ps[p].v

sim.move_to_com()

minmass = 1.e10
for p in ps[1:]:
    if p.m < minmass:
        minmass = p.m
minhill = ps[1].a*(2.*minmass/(3.*ps[0].m))**(1./3.)
for p in ps[1:]:
    p.r = minhill

filename='data/IC{0}K{1:.4e}mag{2:.4e}.bin'.format(simID, K, mag)
sim.initSimulationArchive(filename, interval=1.e3)
sim.integrate(0.)
