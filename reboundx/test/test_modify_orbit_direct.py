import rebound
import reboundx
import unittest
import numpy as np

class TestModifyOrbitDirect(unittest.TestCase):
    def modify_orbit_direct(self, integrator="ias15"):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-6,a=1,e=0.1, inc=0.1)
        sim.add(m=1e-6,a=10,e=0.1, inc=0.1)
        sim.dt = sim.particles[1].P/123.0
        sim.integrator = integrator
        sim.move_to_com() 

        rebx = reboundx.Extras(sim)
        mod = rebx.load_operator("modify_orbits_direct")
        rebx.add_operator(mod)
        
        tmax = 1.e3
        
        sim.particles[1].params["tau_e"] = -tmax/10.
        sim.particles[1].params["tau_inc"] = -tmax/10.

        sim.particles[2].params["tau_e"] = -tmax
        sim.particles[2].params["tau_inc"] = -tmax
   
        sim.integrate(tmax)

        self.assertAlmostEqual(sim.particles[2].e, 0.1*np.exp(tmax/sim.particles[2].params["tau_e"]), places=5)
        self.assertAlmostEqual(sim.particles[1].e, 0.1*np.exp(tmax/sim.particles[1].params["tau_e"]), places=5)
        self.assertAlmostEqual(sim.particles[2].inc, 0.1*np.exp(tmax/sim.particles[2].params["tau_inc"]), places=5)
        self.assertAlmostEqual(sim.particles[1].inc, 0.1*np.exp(tmax/sim.particles[1].params["tau_inc"]), places=5)

    def test_modify_orbit_direct_ias15(self):
        self.modify_orbit_direct("ias15")

    def test_modify_orbit_direct_whfast(self):
        self.modify_orbit_direct("whfast")


if __name__ == '__main__':
    unittest.main()

