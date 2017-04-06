# python run.py 4res taue0 tmax where 4res means the starting binary is in data/4res.bin and taue0 is the initial taue for the innermost planet
import rebound
import reboundx
import numpy as np
import sys

def collision(reb_sim, col):
    reb_sim.contents._status = 5
    return 0

filename = sys.argv[1]
taue0 = float(sys.argv[2])
tmax = float(sys.argv[3])

sim = rebound.Simulation.from_file("data/"+filename+".bin")
ps = sim.particles

sim.collision = "direct"
sim.collision_resolve = collision
sim.ri_whfast.safe_mode = 0

minmass = 1.e10
for p in ps[1:]:
    if p.m < minmass:
        minmass = p.m
minhill = ps[1].a*(2.*minmass/(3.*ps[0].m))**(1./3.)
for p in ps[1:]:
    p.r = minhill

rebx = reboundx.Extras(sim)
rebx.add("tides_synchronous_ecc_damping")

c = 63197.8
gr = rebx.add("gr_potential")
gr.params["c"] = c

aref = ps[1].a
mref = ps[1].m
for p in ps[1:]:
    p.params["tides_synchronous_tau_e"] = -taue0*(p.m/mref)*(p.a/aref)**6.5

rebx.save("data/"+filename+"rebxGR"+"{0:.2e}".format(taue0)+".bin")
sim.initSimulationArchive("data/"+filename+"archiveGR"+"{0:.2e}".format(taue0)+".bin", interval=tmax/1000.)
sim.integrate(tmax)
