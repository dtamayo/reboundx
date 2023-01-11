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
        L = sim.angular_momentum()
        # set up spin L = I Omega to be opposite the orbital angular momentum
        rebx = reboundx.Extras(sim)
        ps[1].params['I'] = 1.
        ps[1].params['Omega'] = [-L[0], -L[1], -L[2]]
        Lspin = rebx.spin_angular_momentum()

        self.assertLess(L[0]+Lspin[0], 1e-15)
        self.assertLess(L[1]+Lspin[1], 1e-15)
        self.assertLess(L[2]+Lspin[2], 1e-15)
    
    def test_rotate_sim(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3, a=1, inc=0.4)
        ps = sim.particles
        L = sim.angular_momentum()
        # set up spin L = I Omega to be opposite the orbital angular momentum
        rebx = reboundx.Extras(sim)
        ps[1].params['I'] = 1.
        ps[1].params['Omega'] = [L[0], L[1], L[2]]
        Lspin = rebx.spin_angular_momentum()
        Ltot = [L[0]+Lspin[0], L[1]+Lspin[1],L[2]+Lspin[2]]
        L = (Ltot[0]**2 + Ltot[1]**2 + Ltot[2]**2)**(0.5)
        rot = rebound.Rotation.to_new_axes(newz=Ltot)
        rebx.rotate_simulation(rot)
        
        Lnew = sim.angular_momentum()
        Lspinnew = rebx.spin_angular_momentum()
        Ltotnew = [Lnew[0]+Lspinnew[0], Lnew[1]+Lspinnew[1],Lnew[2]+Lspinnew[2]]

        self.assertLess(Ltotnew[0], 1e-15)
        self.assertLess(Ltotnew[1], 1e-15)
        self.assertAlmostEqual(Ltotnew[2], L, delta=1e-15)


if __name__ == '__main__':
    unittest.main()
