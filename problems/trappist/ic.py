import rebound
import reboundx
import numpy as np
from collections import OrderedDict
from numpy.random import seed, normal, uniform
import matplotlib.pyplot as plt

def wrap(val):
    while val < 0:
        val += 2*np.pi
    while val > 2*np.pi:
        val -= 2*np.pi
    return val*180/np.pi

def output(sim,planets,resonances,threebodyresonances,outputs):
    t, e, P, Pratio, phi1, phi2, deltapomega, phi3body = outputs
    ps = sim.particles
    for p in planets:
        try:
            e[p].append(ps[p].e)
            P[p].append(ps[p].P)
            t[p].append(sim.t)
        except:
            pass
    for resonance in resonances.items():
        pair = resonance[0]
        p1 = pair[0]
        p2 = pair[1]
        try:
            Pratio[pair].append(ps[p2].P/ps[p1].P)
            t[pair].append(sim.t)
            if resonance[1] is not None:
                res = resonance[1]
                p = res[1]
                q = res[0]-res[1]
                phi1[pair].append(wrap((p+q)*ps[p2].l - p*ps[p1].l - q*ps[p2].pomega))
                phi2[pair].append(wrap((p+q)*ps[p2].l - p*ps[p1].l - q*ps[p1].pomega))
                deltapomega[pair].append(wrap(ps[p2].pomega - ps[p1].pomega))
        except:
            pass
        
    for resonance in threebodyresonances.items():
        triad = resonance[0]
        p1 = triad[0]
        p2 = triad[1]
        p3 = triad[2]
        res = resonance[1]
        p = res[0]
        q = res[1]
        try:
            phi3body[triad].append(wrap(p*ps[p1].l - (p+q)*ps[p2].l + q*ps[p3].l))
            t[triad].append(sim.t)
        except:
            pass

def plot(planets,resonances,threebodyresonances,outputs):
    t, e, P, Pratio, phi1, phi2, deltapomega, phi3body = outputs
    fig, axarr = plt.subplots(ncols=2, nrows=4, figsize=(18,10))
    for p in planets:
        axarr[0,0].plot(t[p], P[p], '.', label=p)
        axarr[0,1].plot(t[p], e[p], '.', label=p)
    for resonance in resonances.items():
        pair = resonance[0]
        res = resonance[1]
        axarr[1,0].plot(t[pair], Pratio[pair], '.', label=pair)
        if res is not None:
            resratio = res[0]/res[1]
            axarr[1,1].plot(t[pair], np.array(Pratio[pair])-resratio, '.', label=pair)
            axarr[2,0].plot(t[pair], phi1[pair], '.', label=pair)
            axarr[2,1].plot(t[pair], phi2[pair], '.', label=pair)
            axarr[3,0].plot(t[pair], deltapomega[pair], '.', label=pair)
            axarr[2,0].set_ylim([0.,360])
            axarr[2,1].set_ylim([0.,360])
    for resonance in threebodyresonances.items():
        triad = resonance[0]
        axarr[3,1].plot(t[triad], phi3body[triad], '.', label=triad)
    
    for ax in axarr.flatten():
        ax.legend(loc='upper right')
    return fig

def removedamping(sim, Tremoval, K, planets, resonances, threebodyresonances, outputs): # Tremovel in # of taue0s
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
        output(sim,planets,resonances,threebodyresonances,outputs)
        
    for p in ps[1:]:
        p.params["tau_e"] = np.inf
        p.params["tau_a"] = np.inf

def integrate(sim, T, planets, resonances, threebodyresonances, outputs):
    Nout = 1000
    T0 = sim.t
    times = np.linspace(T0, T0+T, Nout)
    ps = sim.particles
    for i, time in enumerate(times):
        sim.integrate(time)
        output(sim,planets,resonances,threebodyresonances,outputs)
        
def initialize(planets, resonances, threebodyresonances):
    t, e, P, Pratio, phi1, phi2, deltapomega, phi3body = {}, {}, {}, {}, {}, {}, {}, {}

    for label in planets:
        t[label] = []
        e[label] = []
        P[label] = []

    for pair in resonances.keys():
        t[pair] = []
        Pratio[pair] = []
        phi1[pair] = []
        phi2[pair] = []
        deltapomega[pair] = []
        
    for triad in threebodyresonances.keys():
        t[triad] = []
        phi3body[triad] = []
        
    return (t, e, P, Pratio, phi1, phi2, deltapomega, phi3body)

def drawnormal(pair):
    mean, sigma = pair
    val = -1
    while val < 0:
        val = normal(mean, sigma)
    return val

def res_chain_setup(sim, planets, resonances, delta, massdist, incdist):
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

def plotsa(sa,planets,resonances,threebodyresonances,loc='upper left',istart=0,ifinal=None):
    if ifinal is None:
        ifinal = len(sa)
    Nout = len(sa)
    N = sa[0].N
   
    P0s = []
    for p in sa[0].particles[1:]:
        P0s.append(p.P)
 
    outputs = initialize(planets, resonances, threebodyresonances)
    for i in range(istart, ifinal):
        sim = sa[i]
        ps = sim.particles
        output(sim,planets,resonances,threebodyresonances,outputs)
    
    t, e, P, Pratio, phi1, phi2, deltapomega, phi3body = outputs
        
    fig, axarr = plt.subplots(ncols=2, nrows=4, figsize=(18,10))
    for i, p in enumerate(planets):
        axarr[0,0].plot(t[p], np.array(P[p])-P0s[i], '.', label=p)
        axarr[0,1].plot(t[p], e[p], '.', label=p)
    for resonance in resonances.items():
        pair = resonance[0]
        res = resonance[1]
        resratio = res[0]/res[1]
        axarr[1,0].plot(t[pair], Pratio[pair], '.', label=pair)
        axarr[1,1].plot(t[pair], np.array(Pratio[pair])-resratio, '.', label=pair)
        axarr[2,0].plot(t[pair], phi1[pair], '.', label=pair)
        axarr[2,1].plot(t[pair], phi2[pair], '.', label=pair)
        axarr[3,0].plot(t[pair], deltapomega[pair], '.', label=pair)
        axarr[2,0].set_ylim([0.,360])
        axarr[2,1].set_ylim([0.,360])
    for resonance in threebodyresonances.items():
        triad = resonance[0]
        axarr[3,1].plot(t[triad], phi3body[triad], '.', label=triad)
    for ax in axarr.flatten():
        ax.legend(loc=loc)
        
    return fig
