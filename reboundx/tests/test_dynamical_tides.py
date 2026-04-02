import rebound
import reboundx
import unittest
import os
import numpy as np

class TestTides(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.G = 4*np.pi**2
        self.sim.add(m=1.)
        self.sim.add(m=1.e-3, a=1.5, e=0.987, r=7.5e-4)
        self.sim.move_to_com()
        self.rebx = reboundx.Extras(self.sim)
        self.force = self.rebx.load_force("tides_dynamical")
        self.rebx.add_force(self.force)

    def test_noM(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1., r=0.005)
        self.sim.add(a=0.2, e=0.1, r=1e-4)
        self.sim.move_to_com()
        self.rebx = reboundx.Extras(self.sim)
        self.force = self.rebx.load_force("tides_dynamical")
        self.rebx.add_force(self.force)
        with self.assertRaises(RuntimeError):
            self.sim.integrate(10)

    def test_noR(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1., r=0.005)
        self.sim.add(m=1., a=0.2, e=0.1)
        self.sim.move_to_com()
        self.rebx = reboundx.Extras(self.sim)
        self.force = self.rebx.load_force("tides_dynamical")
        self.rebx.add_force(self.force)
        with self.assertRaises(RuntimeError):
            self.sim.integrate(10)

    def test_a_decay(self):
        self.sim.integrate(1e3)
        self.assertLess(self.sim.particles[1].a, 1.)

    def test_e_decay(self):
        self.sim.integrate(1e3)
        self.assertLess(self.sim.particles[1].e, 0.98)

if __name__ == '__main__':
    unittest.main()
