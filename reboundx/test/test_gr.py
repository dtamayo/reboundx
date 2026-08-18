import rebound
import reboundx
import unittest
import os
import numpy as np
import reboundx.constants as constants

class TestGR(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1., r=0.005)
        self.sim.add(m=1.e-3, a=1.0, e=0.1)
        self.sim.dt = self.sim.particles[1].P/20.2
        self.sim.move_to_com()
        self.rebx = reboundx.Extras(self.sim)

    def test_gr_potential_error(self):
        self.force = self.rebx.load_force("gr_potential")
        self.rebx.add_force(self.force)
        with self.assertRaises(RuntimeError):
            self.sim.step() # didn't set c

    def test_gr_error(self):
        self.force = self.rebx.load_force("gr")
        self.rebx.add_force(self.force)
        with self.assertRaises(RuntimeError):
            self.sim.step() # didn't set c

    def test_gr_full_error(self):
        self.force = self.rebx.load_force("gr_full")
        self.rebx.add_force(self.force)
        with self.assertRaises(RuntimeError):
            self.sim.step() # didn't set c

    def test_gr_potential_energy(self):
        self.force = self.rebx.load_force("gr_potential")
        self.rebx.add_force(self.force)
        self.force.params['c'] = constants.C
        E0 = self.sim.energy() + self.rebx.gr_potential_potential(self.force)
        self.sim.integrate(1e4)
        E = self.sim.energy() + self.rebx.gr_potential_potential(self.force)
        self.assertLess(abs((E-E0)/E0), 1e-12)

    def test_gr_energy(self):
        self.force = self.rebx.load_force("gr")
        self.rebx.add_force(self.force)
        self.force.params['c'] = constants.C
        E0 = self.rebx.gr_hamiltonian(self.force)
        self.sim.integrate(1e4)
        E = self.rebx.gr_hamiltonian(self.force)
        self.assertLess(abs((E-E0)/E0), 1e-12)

    def test_gr_full_energy(self):
        self.force = self.rebx.load_force("gr_full")
        self.rebx.add_force(self.force)
        self.force.params['c'] = constants.C
        E0 = self.rebx.gr_full_hamiltonian(self.force)
        self.sim.integrate(1e4)
        E = self.rebx.gr_full_hamiltonian(self.force)
        self.assertLess(abs((E-E0)/E0), 1e-12)

    def test_Nactive(self):
        sim = rebound.Simulation()
        sim.add(m=1., r=0.005)
        sim.add(m=1.e-3, a=1.0, e=0.1)
        sim.add(a=2.0, e=0.1)
        sim.add(a=3.0, e=0.1)
        sim.dt = sim.particles[1].P/20.2
        sim.move_to_com()
        rebx = reboundx.Extras(sim)
        force = rebx.load_force("gr")
        rebx.add_force(force)
        force.params['c'] = constants.C
        sim.integrate(100.)

        sim2= rebound.Simulation()
        sim2.add(m=1., r=0.005)
        sim2.add(m=1.e-3, a=1.0, e=0.1)
        sim2.add(a=2.0, e=0.1)
        sim2.add(a=3.0, e=0.1)
        sim2.dt = sim2.particles[1].P/20.2
        sim2.move_to_com()
        rebx2 = reboundx.Extras(sim2)
        force2 = rebx2.load_force("gr")
        rebx2.add_force(force2)
        force2.params['c'] = constants.C
        sim2.N_active=2
        sim2.integrate(100.)

        x0 = sim.particles[3].x
        x = sim2.particles[3].x
        self.assertLess(abs((x-x0)/x0), 1e-12)

    # add energy errors and gr gives same output with and without Nactive
if __name__ == '__main__':
    unittest.main()
