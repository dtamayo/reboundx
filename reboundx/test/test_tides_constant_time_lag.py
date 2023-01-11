import rebound
import reboundx
import unittest
import os
import numpy as np

class TestTides(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1., r=0.005)
        self.sim.add(m=1.e-3, a=0.05, e=0.1, r=0.0005)
        self.sim.move_to_com()
        self.rebx = reboundx.Extras(self.sim)
        self.force = self.rebx.load_force("tides_constant_time_lag")
        self.rebx.add_force(self.force)

    def test_conservation_star_highmratio(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1., r=0.005)
        self.sim.add(m=1., a=0.2, e=0.1, r=0.005)
        self.sim.move_to_com()
        self.rebx = reboundx.Extras(self.sim)
        self.force = self.rebx.load_force("tides_constant_time_lag")
        self.rebx.add_force(self.force)
        ps = self.sim.particles
        ps[0].params['tctl_k2'] = 0.04
        self.do_test_conservation()

    def test_conservation_planet_highmratio(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1., r=0.005)
        self.sim.add(m=1., a=0.2, e=0.1, r=0.005)
        self.sim.move_to_com()
        self.rebx = reboundx.Extras(self.sim)
        self.force = self.rebx.load_force("tides_constant_time_lag")
        self.rebx.add_force(self.force)
        ps = self.sim.particles
        ps[1].params['tctl_k2'] = 0.04
        self.do_test_conservation()
    
    def test_conservation_star(self):
        ps = self.sim.particles
        ps[0].params['tctl_k2'] = 0.04
        self.do_test_conservation()

    def test_conservation_planet(self):
        ps = self.sim.particles
        ps[1].params['tctl_k2'] = 0.4
        self.do_test_conservation()
    
    def test_conservation_star_movecom(self):
        # if you go much larger, IAS15 starts giving errors due to roundoff
        self.sim.particles[0].vy += 0.01 
        self.sim.particles[1].vy += 0.01 
        self.test_conservation_star()
    
    def test_conservation_planet_movecom(self):
        self.sim.particles[0].vy += 0.01 
        self.sim.particles[1].vy += 0.01 
        self.test_conservation_planet()

    def do_test_conservation(self):
        H0 = self.sim.energy() + self.rebx.tides_constant_time_lag_potential(self.force)
        self.sim.integrate(1.e3*self.sim.particles[1].P)
        H = self.sim.energy() + self.rebx.tides_constant_time_lag_potential(self.force)
        self.assertLess(abs((H-H0)/H0), 1.e-12)
        self.assertGreater(self.sim.particles[1].pomega, 1.e-6)

class TestTidesAnalytic(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=0.86, r = 0.78)   
        self.sim.add(m=3.e-6, a=1., e=0.05)
        self.sim.move_to_com()
        ps = self.sim.particles

        self.rebx = reboundx.Extras(self.sim)
        self.tides = self.rebx.load_force("tides_constant_time_lag")
        self.rebx.add_force(self.tides)
        ps[0].params["tctl_k2"] = 0.023 # in AU
        ps[0].params["tctl_tau"] = 0.3
        ps[0].params["Omega"] = 0.

        self.q = (ps[1].m/ps[0].m)
        self.T = ps[0].r**3/self.sim.G/ps[0].m/ps[0].params["tctl_tau"]
        self.taua = self.T/6/ps[0].params["tctl_k2"]/self.q/(1+self.q)*(ps[1].a/ps[0].r)**8

    def test_adamping(self):
        ps = self.sim.particles
        tmax = 4e4*ps[1].P
        apred = ps[0].r*((ps[1].a/ps[0].r)**8 - 48.*ps[0].params["tctl_k2"]*self.q*(1+self.q)*tmax/self.T)**(1./8.)

        self.sim.integrate(tmax)
        self.assertLess(abs((ps[1].a-apred)/apred), 1.e-2) # 1%
    
    def test_linear_adamping(self):
        ps = self.sim.particles
        tmax = self.taua/100
        apred = ps[1].a*np.exp(-tmax/self.taua)

        self.sim.integrate(tmax)
        self.assertLess(abs((ps[1].a-apred)/apred), 0.1) # 10%
    
    def test_linear_edamping(self):
        ps = self.sim.particles
        tmax = self.taua/100
        epred = ps[1].e*np.exp(-tmax/(6./27.*self.taua))

        self.sim.integrate(tmax)
        self.assertLess(abs((ps[1].e-epred)/epred), 0.1) # 10%
    
    def test_adamping_movecom(self):
        self.sim.particles[0].vy += 0.01 
        self.sim.particles[1].vy += 0.01 
        self.test_adamping()
    
    def test_linear_adamping_movecom(self):
        self.sim.particles[0].vy += 0.01 
        self.sim.particles[1].vy += 0.01 
        self.test_linear_adamping()
    
    def test_linear_edamping_movecom(self):
        self.sim.particles[0].vy += 0.01 
        self.sim.particles[1].vy += 0.01 
        self.test_linear_edamping()

if __name__ == '__main__':
    unittest.main()
