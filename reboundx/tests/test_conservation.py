import rebound
import reboundx
import unittest
import os

class TestConservation(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(m=1.e-4, a=1., e=0.1)
        self.sim.add(m=1.e-4, a=2, e=0.1, inc=0.2)
        self.sim.integrator="whfast"
        self.sim.dt = self.sim.particles[1].P/100
        self.sim.move_to_com()
    
    def test_gr_full(self):
        name = 'gr_full'
        sim = self.sim
        sim.integrator = "ias15"
        rebx = reboundx.Extras(sim)
        force = rebx.load_force(name)
        rebx.add_force(force)
        force.params['c'] = 1.e4
        ps = sim.particles
        H0 = rebx.gr_full_hamiltonian(force)
        sim.integrate(1.e4)
        H = rebx.gr_full_hamiltonian(force)
        self.assertLess(abs((H-H0)/H0), 1.e-12)
    
    def test_gr(self):
        name = 'gr'
        sim = self.sim
        sim.integrator = "ias15"
        rebx = reboundx.Extras(sim)
        force = rebx.load_force(name)
        rebx.add_force(force)
        force.params['c'] = 1.e4
        ps = sim.particles
        H0 = rebx.gr_hamiltonian(force)
        sim.integrate(1.e4)
        H = rebx.gr_hamiltonian(force)
        self.assertLess(abs((H-H0)/H0), 1.e-12)
    
    def test_gr_potential(self):
        name = 'gr_potential'
        sim = self.sim
        sim.integrator = "ias15"
        rebx = reboundx.Extras(sim)
        force = rebx.load_force(name)
        rebx.add_force(force)
        force.params['c'] = 1.e4
        ps = sim.particles
        H0 = sim.energy() + rebx.gr_potential_potential(force)
        sim.integrate(1.e4)
        H = sim.energy() + rebx.gr_potential_potential(force)
        self.assertLess(abs((H-H0)/H0), 1.e-12)

    def test_central_force(self):
        name = 'central_force'
        sim = self.sim
        sim.integrator = "ias15"
        rebx = reboundx.Extras(sim)
        force = rebx.load_force(name)
        rebx.add_force(force)
        ps = sim.particles
        ps[0].params['Acentral'] = 1.e-4
        ps[0].params['gammacentral'] = -1
        H0 = sim.energy() + rebx.central_force_potential()
        sim.integrate(1.e4)
        H = sim.energy() + rebx.central_force_potential()
        self.assertLess(abs((H-H0)/H0), 1.e-12)

    def test_gravitational_harmonics(self):
        name = 'gravitational_harmonics'
        sim = self.sim
        sim.integrator = "ias15"
        rebx = reboundx.Extras(sim)
        force = rebx.load_force(name)
        rebx.add_force(force)
        ps = sim.particles
        ps[0].params['J2'] = 1.e-3
        ps[0].params['J4'] = 1.e-3
        ps[0].params['R_eq'] = 1.e-3
        H0 = sim.energy() + rebx.gravitational_harmonics_potential()
        sim.integrate(1.e4)
        H = sim.energy() + rebx.gravitational_harmonics_potential()
        self.assertLess(abs((H-H0)/H0), 1.e-12)

if __name__ == '__main__':
    unittest.main()
