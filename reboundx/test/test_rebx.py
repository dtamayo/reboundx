import rebound
import reboundx
import unittest

class TestRebx(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(a=1., e=0.2)
        self.rebx = reboundx.Extras(self.sim)

    def test_simple(self):
        self.gr = self.rebx.load_force('gr')
        self.rebx.add_force(self.gr)
        self.gr.params['c'] = 1e2

        self.sim.integrate(10)
        self.assertGreater(self.sim.particles[1].pomega, 1.e-4)

if __name__ == '__main__':
    unittest.main()
