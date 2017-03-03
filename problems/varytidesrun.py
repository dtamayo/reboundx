import numpy as np
import rebound
import reboundx
import sys
from random import random, uniform, seed

P = [1.51087081, 2.4218233, 4.049610, 6.099615, 9.206690, 12.35294, 20] # days
P = [i/365.25 for i in P] # years
a = [11.11, 15.21, 21.44, 28.17, 37.1, 45.1, 63]
a = [i*1.e-3 for i in a] # AU
scalepar = [20.50, 28.08, 39.55, 51.97, 68.4, 83.2, 117] # dimensionless a/Rstar
radius = [1.086, 1.056, 0.772, 0.918, 1.045, 1.127, 0.755] # Rearth
RearthinAU = 4.26e-5
m = [0.85, 1.38, 0.41, 0.62, 0.68, 1.34] # Mearth
earthmass = 3.e-6 # solar masses
rho = [0.66, 1.17, 0.89, 0.8, 0.6, 0.94] # rhoEarth
radiustostar = [0.7266, 0.687, 0.367, 0.519, 0.673, 0.782, 0.352]
radiustostar = [np.sqrt(i/100.) for i in radiustostar]
Q = 100
tsync = [20*Q*rho[i]*scalepar[i]**3*P[i]**2 for i in range(6)]
tau_e = [1.5e-6*Q*scalepar[i]**5/radiustostar[i]**5*m[i]*P[i] for i in range(6)]
tau_a = [i*1e4 for i in tau_e]

G = 6.67e-11 # all SI
Mstar = 0.08*2e30 
earthinkg = 6e24 
earthradius = 6400000
AUinm = 1.5e11
yrins = 3.15e7
earthheat = 4.7e13
solarconstant = 1361
Edot = [G*Mstar*m[i]*earthinkg/2/(a[i]*AUinm)/(tau_a[i]*yrins)/earthheat for i in range(6)] # compared to Earth heating rate
Stidal = [Edot[i]*earthheat/np.pi/(radius[i]*earthradius)**2/solarconstant for i in range(6)] # compared to solar constant
irradiation = [4.25, 2.27, 1.143, 0.662, 0.382, 0.258, 0.131]
internaltoirradiation = [Stidal[i]/irradiation[i] for i in range(6)]
hillradius = [a[i]*(earthinkg*m[i]/(3.*Mstar))**(1/3) for i in range(6)]

def collision(reb_sim, col):
    reb_sim.contents._status = 5
    #for p in reb_sim.contents.particles[1:]:
    #    print(p.orbit)
    return 0

sim = rebound.Simulation()
sim.G = 4*np.pi**2

e0 = 0.01
i0 = 0.5*np.pi/180.
sim_id = sys.argv[1]
tauefac = float(sys.argv[2])
seed(sim_id)

rebx = reboundx.Extras(sim)
edamp = rebx.add("tides_synchronous_ecc_damping")

sim.add(m=0.08)
ps = sim.particles

for i in range(6):
    sim.add(m=m[i]*earthmass,P=P[i], e=e0, pomega=random()*2.*np.pi, inc=i0, Omega=random()*2.*np.pi, f=random()*2.*np.pi, r=hillradius[i])
    ps[i+1].params["tides_synchronous_tau_e"] = tau_e[i]*tauefac
sim.move_to_com() # Moves to the center of momentum frame

sim.integrator="whfast"
sim.dt=sim.particles[1].P*0.09
sim.collision = "direct"
sim.collision_resolve = collision

tmax = 1.e5*ps[1].P

sim.integrate(tmax)

with open("tidesdata/run"+str(sim_id)+".txt", "w") as f:
    f.write("{0},{1}".format(sim.t, tauefac))
    for i in range(1,sim.N):
        f.write(",{0}".format(ps[i].P))
    f.write("\n")
