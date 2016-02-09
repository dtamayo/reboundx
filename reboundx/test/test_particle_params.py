import rebound
import reboundx
from reboundx import data
import unittest
import math

class Test_particle_params(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        data.add_earths(self.sim, ei=1.e-3) 
        self.TINY=1.e-15 
        rebx = reboundx.Extras(self.sim)

    def tearDown(self):
        self.sim = None
    
    def test_particle_params(self):
        """
        Assign all the particle to two particles, and make sure we get the same values back
        """
        p1 = self.sim.particles[1]
        p2 = self.sim.particles[2]
        param = "tau_a"

        if hasattr(p1, param): # make sure we're not just adding a new member to the particle and getting it back
            setattr(p1, param, float(3))
        if hasattr(p2, param):
            setattr(p2, param, float(5))
        self.assertAlmostEqual(getattr(p1, param), 3, delta=self.TINY)
        self.assertAlmostEqual(getattr(p2, param), 5, delta=self.TINY)

if __name__ == '__main__':
    unittest.main()
