import sys

import rebound
import reboundx
import numpy as np

ID = sys.argv[1]

Nout = 1000
maxorb = 1.e3
C2 = 1.e6
dt = 7.e-2

def run(name, rebintegrator, rebxintegrator, order, force_as_operator, ID):
    sim = rebound.Simulation()
    sim.add(m=1.)
    sim.add(m=1.e-3, a=1., e=0.1)
    #sim.add(m=1.e-3, a=2., e=0.1)
    sim.move_to_com()
    sim.integrator=rebintegrator
    sim.dt = dt*sim.particles[1].P
    rebx = reboundx.Extras(sim)
    rebx.integrator = rebxintegrator
    gr = rebx.add("gr")
    gr.params["c"] = np.sqrt(C2)
    gr.operator_order = order
    gr.force_as_operator = force_as_operator
    Eerr = np.zeros(Nout)
    E0 = rebx.gr_hamiltonian(sim, gr)
    times = np.logspace(0,np.log10(maxorb*sim.particles[1].P),Nout)

    for i, time in enumerate(times):
        sim.integrate(time, exact_finish_time=0)
        E = rebx.gr_hamiltonian(sim, gr)
        with open("longtermdata/"+name+str(ID)+".txt", "w") as f:
            f.write("{0:e}\t{1:e}\t{2:.16f}\t{3:.16f}\n".format(sim.t, abs((E-E0)/E0), sim.particles[1].a, sim.particles[1].e))

run("euler", "whfast", "euler", 2, 1, ID)
run("implicit_midpoint", "whfast", "implicit_midpoint", 2, 1, ID)
run("none", "whfast", "none", 2, 1, ID)
run("naive", "whfast", "none", 2, 0, ID)
