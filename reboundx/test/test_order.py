import rebound
import reboundx
import unittest
import numpy as np

class TestOrder(unittest.TestCase):
    def setUp(self):
        Mp = 1.e-8
        Nout = 1000
        self.J2 = 1.e-4
        self.dtOverT = 0.01
        self.times = np.logspace(0, 4, 1000)

        self.sim = rebound.Simulation()
        self.sim.G = 4*np.pi**2
        self.sim.add(m=1.)
        self.sim.add(m=Mp, P=1., e=0.1)
        self.sim.add(m=Mp, P=2.3, e=0.1)
        self.sim.move_to_com()
        self.sim.integrator="whfast"
        self.sim.dt = self.dtOverT*self.sim.particles[1].P

        self.rebx = reboundx.Extras(self.sim)
        self.rebx.integrator = 'rk4'
        self.gh = self.rebx.add("gravitational_harmonics")
        ps = self.sim.particles
        ps[0].params['J2'] = self.J2
        ps[0].params['R_eq'] = ps[1].a
        self.gh.force_as_operator = 1

    def tearDown(self):
        self.sim = None
        self.rebx = None
        self.gh = None
        self.J2 = None
        self.dtOverT = None
        self.times = None
 
    def test_1storder(self):
        self.gh.operator_order = 1
        Eerr = np.zeros(1000)
        E0 = self.sim.calculate_energy() + self.rebx.gravitational_harmonics_hamiltonian(self.sim)

        for i, time in enumerate(self.times):
            self.sim.integrate(time, exact_finish_time=0)
            E = self.sim.calculate_energy() + self.rebx.gravitational_harmonics_hamiltonian(self.sim)
            Eerr[i] = np.abs((E-E0)/E0)
        
        self.assertLess(np.max(Eerr), np.pi*self.J2*self.dtOverT)

    def test_2ndorder(self):
        self.gh.operator_order = 2
        Eerr = np.zeros(1000)
        E0 = self.sim.calculate_energy() + self.rebx.gravitational_harmonics_hamiltonian(self.sim)

        for i, time in enumerate(self.times):
            self.sim.integrate(time, exact_finish_time=0)
            E = self.sim.calculate_energy() + self.rebx.gravitational_harmonics_hamiltonian(self.sim)
            Eerr[i] = np.abs((E-E0)/E0)
        
        self.assertLess(np.max(Eerr), np.pi**2/6.*self.J2*self.dtOverT**2)

if __name__ == '__main__':
    unittest.main()
