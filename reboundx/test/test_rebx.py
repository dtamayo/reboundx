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
    
    def test_total_angular_momentum(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3, a=1, inc=0.4)
        ps = sim.particles
        L = sim.calculate_angular_momentum()
        # set up spin L = Iw = moi*s to be opposite the orbital angular momentum
        rebx = reboundx.Extras(sim)
        ps[1].params['moi'] = 1.
        ps[1].params['sx'] = -L[0]
        ps[1].params['sy'] = -L[1]
        ps[1].params['sz'] = -L[2]
        Ltot = rebx.calculate_total_angular_momentum()
        self.assertLess(Ltot[0], 1e-15)
        self.assertLess(Ltot[1], 1e-15)
        self.assertLess(Ltot[2], 1e-15)
    
    def test_rotate_sim(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3, a=1, inc=0.4)
        ps = sim.particles
        L = sim.calculate_angular_momentum()
        # set up spin L = Iw = moi*s to be opposite the orbital angular momentum
        rebx = reboundx.Extras(sim)
        ps[1].params['moi'] = 1.
        ps[1].params['sx'] = L[0]
        ps[1].params['sy'] = L[1]
        ps[1].params['sz'] = L[2]
        Ltot = rebx.calculate_total_angular_momentum()
        L = (Ltot[0]**2 + Ltot[1]**2 + Ltot[2]**2)**(0.5)
        rot = rebound.Rotation.to_new_axes(newz=Ltot)
        rebx.rotate_simulation(rot)
        Ltotnew = rebx.calculate_total_angular_momentum()
        self.assertLess(Ltotnew[0], 1e-15)
        self.assertLess(Ltotnew[1], 1e-15)
        self.assertAlmostEqual(Ltotnew[2], L, delta=1e-15)


if __name__ == '__main__':
    unittest.main()
