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

    def test_detach(self):
        self.gr = self.rebx.load_force('gr')
        self.rebx.add_force(self.gr)
        self.gr.params['c'] = 1e2

        self.rebx.detach()
        self.sim.integrate(10)
        self.assertLess(self.sim.particles[1].pomega, 1.e-10) # should be 0

    def test_create_in_function(self):
        def makesim():
            sim = rebound.Simulation()
            sim.add(m=1.)
            sim.add(a=1.)
            rebx = reboundx.Extras(sim)
            gr = rebx.load_force('gr')
            rebx.add_force(gr)
            gr.params['c'] = 1e2
            return sim, rebx

        sim, rebx = makesim()
        sim.integrate(10)
        self.assertGreater(sim.particles[1].pomega, 0.01)

    def test_returnsim(self):
        def returnsim():
            sim = rebound.Simulation()
            sim.add(m=1.)
            sim.add(a=1.)
            rebx = reboundx.Extras(sim)
            gr = rebx.load_force('gr')
            rebx.add_force(gr)
            gr.params['c'] = 1e2
            return sim

        sim = returnsim()
        sim.integrate(10)
        self.assertGreater(sim.particles[1].pomega, 0.01)

    def test_delete(self):
        def makesim():
            sim = rebound.Simulation()
            sim.add(m=1.)
            sim.add(a=1.)
            rebx = reboundx.Extras(sim)
            gr = rebx.load_force('gr')
            rebx.add_force(gr)
            gr.params['c'] = 1e2
            return sim, rebx

        for i in range(10):
            sim, rebx = makesim()
            sim.integrate(10)

        self.assertGreater(sim.particles[1].pomega, 0.01)

    def test_rebx_not_attached(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["a"] = 7

if __name__ == '__main__':
    unittest.main()
