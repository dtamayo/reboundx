import rebound
import reboundx
import unittest
import os
import numpy as np

class TestStochastic(unittest.TestCase):
    def test_default(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1., r=0.005)
        self.sim.add(m=1.e-3, a=1.0, e=0.1)
        self.sim.integrator = "whfast"
        self.sim.rand_seed = 0
        self.sim.dt = self.sim.particles[1].P/20.2
        self.sim.move_to_com()
        self.rebx = reboundx.Extras(self.sim)
        self.force = self.rebx.load_force("stochastic_forces")
        self.rebx.add_force(self.force)
        ps = self.sim.particles
        ps[1].params['kappa'] = 1e-5
        self.sim.integrate(ps[1].P*1000)
        self.assertLess(np.abs(1.-ps[1].a), 0.004)

    def test_scaling(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1., r=0.005)
        self.sim.add(m=1.e-3, a=0.001, e=0.1)
        self.sim.integrator = "whfast"
        self.sim.rand_seed = 0
        self.sim.dt = self.sim.particles[1].P/20.2
        self.sim.move_to_com()
        self.rebx = reboundx.Extras(self.sim)
        self.force = self.rebx.load_force("stochastic_forces")
        self.rebx.add_force(self.force)
        ps = self.sim.particles
        ps[1].params['kappa'] = 1e-5
        self.sim.integrate(ps[1].P*1000)
        self.assertLess(np.abs(0.001-ps[1].a), 0.000004)

    def test_binary(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1., r=0.005)
        self.sim.add(m=1.e-3, a=1, e=0.1)
        self.sim.integrator = "whfast"
        self.sim.rand_seed = 0
        self.sim.dt = self.sim.particles[1].P/20.2
        self.sim.move_to_com()
        self.rebx = reboundx.Extras(self.sim)
        self.force = self.rebx.load_force("stochastic_forces")
        self.rebx.add_force(self.force)
        ps = self.sim.particles
        ps[1].params['kappa'] = 1e-5
        self.sim.integrate(ps[1].P*10)

        self.sim.save("binary.bin")
        self.rebx.save("binary.rebx")

        self.sim.integrate(ps[1].P*20.5)

        sim2 = rebound.Simulation("binary.bin")
        rebx2 = reboundx.Extras(sim2, "binary.rebx")
        
        sim2.integrate(self.sim.t)
        self.assertEqual(self.sim.particles[1].x, sim2.particles[1].x)
        self.assertEqual(self.sim.particles[1].params["kappa"], sim2.particles[1].params["kappa"])

if __name__ == '__main__':
    unittest.main()
