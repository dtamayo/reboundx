import rebound
import reboundx
import unittest
import numpy as np

"""
Acts as both a test on various integration options for forces working (add_force vs step before/after/both)
and for reproducibility with the simulationarchive under those various modes.
"""

# When the scheme error with leapfrog is larger than the perturbation, we get qualitatively wrong results, so don't test
# Cna only add effects as forces to mercurius, so only test individually with add_force
integrators = ['ias15', 'whfast']
rebxintegrators = ['euler', 'rk2', 'rk4', 'implicit_midpoint']

class TestSimulationArchive(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(m=1.e-4, a=1., e=0.2)
        self.sim.move_to_com()
        self.sim.dt = 1.e-2*self.sim.particles[1].P
        self.rebx = reboundx.Extras(self.sim)
        self.gr = self.rebx.load_force('gr')
        self.gr.params['c'] = 1e2
        self.pomega0 = self.sim.particles[1].pomega
        self.E0 = self.rebx.gr_hamiltonian(self.gr)

    def test_force(self):
        for integrator in integrators + ['mercurius']:
            for rebxintegrator in rebxintegrators:
                self.setUp()
                self.sim.integrator = integrator
                self.rebx.add_force(self.gr)
            
                self.sim.integrate(1000)
                E = self.rebx.gr_hamiltonian(self.gr)
                # Test effect is actually doing something to the pericenter
                self.assertGreater(self.sim.particles[1].pomega, 0.1, msg='REB integrator: {0}, REBX integrator: {1}'.format(integrator, rebxintegrator))
                # test energy conservtaion
                self.assertLess(np.abs((E-self.E0)/self.E0), 1.e-4, msg='REB integrator: {0}, REBX integrator: {1}'.format(integrator, rebxintegrator))
               
                # test bitwise reproducibility starting from an intermediate snapshot
                self.sim.simulationarchive_snapshot('test.sa', deletefile=True)
                self.rebx.save('test.rebx')
                self.sim.integrate(2000)
                self.sim.simulationarchive_snapshot('test.sa')

                sa = reboundx.SimulationArchive('test.sa', 'test.rebx')
                simf, rebxf = sa[-1]
                tmax = simf.t
                sim, rebx = sa[0]
                sim.integrate(tmax)
                self.assertEqual(self.sim.particles[1].x, sim.particles[1].x, msg='REB integrator: {0}, REBX integrator: {1}'.format(integrator, rebxintegrator))

    def test_pre(self):
        for integrator in integrators:
            for rebxintegrator in rebxintegrators:
                self.setUp()
                self.sim.integrator = integrator
                self.sim.ri_ias15.epsilon = 0 # use fixed timesteps for ias15 if used, no problem otherwise
                self.integforce = self.rebx.load_operator("integrate_force")
                self.integforce.params['force'] = self.gr
                self.rebx.add_operator(self.integforce, dtfraction=1., timing="pre")
            
                self.sim.integrate(1000)
                E = self.rebx.gr_hamiltonian(self.gr)
                # Test effect is actually doing something to the pericenter
                self.assertGreater(self.sim.particles[1].pomega, 0.1, msg='REB integrator: {0}, REBX integrator: {1}'.format(integrator, rebxintegrator))
                # test energy conservtaion
                self.assertLess(np.abs((E-self.E0)/self.E0), 1.e-4, msg='REB integrator: {0}, REBX integrator: {1}'.format(integrator, rebxintegrator))
               
                # test bitwise reproducibility starting from an intermediate snapshot
                self.sim.simulationarchive_snapshot('test.sa', deletefile=True)
                self.rebx.save('test.rebx')
                self.sim.integrate(2000)
                self.sim.simulationarchive_snapshot('test.sa')

                sa = reboundx.SimulationArchive('test.sa', 'test.rebx')
                simf, rebxf = sa[-1]
                tmax = simf.t
                sim, rebx = sa[0]
                sim.integrate(tmax)
                self.assertEqual(self.sim.particles[1].x, sim.particles[1].x, msg='REB integrator: {0}, REBX integrator: {1}'.format(integrator, rebxintegrator))

    def test_post(self):
        for integrator in integrators:
            for rebxintegrator in rebxintegrators:
                self.setUp()
                self.sim.integrator = integrator
                self.sim.ri_ias15.epsilon = 0 # use fixed timesteps for ias15 if used, no problem otherwise
                self.integforce = self.rebx.load_operator("integrate_force")
                self.integforce.params['force'] = self.gr
                self.rebx.add_operator(self.integforce, dtfraction=1., timing="post")
            
                self.sim.integrate(1000)
                E = self.rebx.gr_hamiltonian(self.gr)
                # Test effect is actually doing something to the pericenter
                self.assertGreater(self.sim.particles[1].pomega, 0.1, msg='REB integrator: {0}, REBX integrator: {1}'.format(integrator, rebxintegrator))
                # test energy conservtaion
                self.assertLess(np.abs((E-self.E0)/self.E0), 1.e-4, msg='REB integrator: {0}, REBX integrator: {1}'.format(integrator, rebxintegrator))
               
                # test bitwise reproducibility starting from an intermediate snapshot
                self.sim.simulationarchive_snapshot('test.sa', deletefile=True)
                self.rebx.save('test.rebx')
                self.sim.integrate(2000)
                self.sim.simulationarchive_snapshot('test.sa')

                sa = reboundx.SimulationArchive('test.sa', 'test.rebx')
                simf, rebxf = sa[-1]
                tmax = simf.t
                sim, rebx = sa[0]
                sim.integrate(tmax)
                self.assertEqual(self.sim.particles[1].x, sim.particles[1].x, msg='REB integrator: {0}, REBX integrator: {1}'.format(integrator, rebxintegrator))

    def test_2ndorder(self):
        for integrator in integrators:
            for rebxintegrator in rebxintegrators:
                self.setUp()
                self.sim.integrator = integrator
                self.sim.ri_ias15.epsilon = 0 # use fixed timesteps for ias15 if used, no problem otherwise
                self.integforce = self.rebx.load_operator("integrate_force")
                self.integforce.params['force'] = self.gr
                self.rebx.add_operator(self.integforce)
            
                self.sim.integrate(1000)
                E = self.rebx.gr_hamiltonian(self.gr)
                # Test effect is actually doing something to the pericenter
                self.assertGreater(self.sim.particles[1].pomega, 0.1, msg='REB integrator: {0}, REBX integrator: {1}'.format(integrator, rebxintegrator))
                # test energy conservtaion
                self.assertLess(np.abs((E-self.E0)/self.E0), 1.e-4, msg='REB integrator: {0}, REBX integrator: {1}'.format(integrator, rebxintegrator))
               
                # test bitwise reproducibility starting from an intermediate snapshot
                self.sim.simulationarchive_snapshot('test.sa', deletefile=True)
                self.rebx.save('test.rebx')
                self.sim.integrate(2000)
                self.sim.simulationarchive_snapshot('test.sa')

                sa = reboundx.SimulationArchive('test.sa', 'test.rebx')
                simf, rebxf = sa[-1]
                tmax = simf.t
                sim, rebx = sa[0]
                sim.integrate(tmax)
                self.assertEqual(self.sim.particles[1].x, sim.particles[1].x, msg='REB integrator: {0}, REBX integrator: {1}'.format(integrator, rebxintegrator))

if __name__ == '__main__':
    unittest.main()

