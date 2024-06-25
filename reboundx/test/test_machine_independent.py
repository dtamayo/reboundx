import rebound
import reboundx
import unittest
import os

class TestMachineIndependent(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(a=1., e=0.2)
        self.rebx = reboundx.Extras(self.sim)

        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1.e-4, a=1., e=0.1)
        sim.add(m=1.e-4, a=2, e=0.1, inc=0.2)
        sim.integrator="whfast"
        sim.dt = sim.particles[1].P/100
        sim.move_to_com()
        sim.save_to_file('binaries/twoplanets.bin', delete_file=True)

    def test_modify_orbits_forces(self):
        sim = rebound.Simulation("binaries/twoplanets.bin")
        sim.save_to_file('binaries/modify_orbits_forces.sa', interval=1e3, delete_file=True)
        rebx = reboundx.Extras(sim)
        mod = rebx.load_force('modify_orbits_forces')
        rebx.add_force(mod)
        ps = sim.particles
        ps[1].params['tau_a'] = -1e4
        ps[1].params['tau_e'] = -1e3
        ps[2].params['tau_e'] = -1e3
        ps[2].params['tau_inc'] = -1e3
        rebx.save('binaries/modify_orbits_forces.rebx')
        sim.integrate(1.e4)
        
        sa = reboundx.Simulationarchive('binaries/modify_orbits_forces.sa', 'binaries/modify_orbits_forces.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)

    def test_gr_full(self):
        sim = rebound.Simulation('binaries/twoplanets.bin')
        sim.save_to_file('binaries/gr_full.sa', interval=1e3, delete_file=True)
        rebx = reboundx.Extras(sim)
        force = rebx.load_force('gr_full')
        rebx.add_force(force)
        force.params['c'] = 1.e4
        ps = sim.particles
        rebx.save('binaries/gr_full.rebx')
        sim.integrate(1.e4)
        sa = reboundx.Simulationarchive('binaries/gr_full.sa', 'binaries/gr_full.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)
    
    def test_gr(self):
        sim = rebound.Simulation('binaries/twoplanets.bin')
        sim.save_to_file('binaries/gr.sa', interval=1e3, delete_file=True)
        rebx = reboundx.Extras(sim)
        force = rebx.load_force('gr')
        rebx.add_force(force)
        force.params['c'] = 1.e4
        ps = sim.particles
        rebx.save('binaries/gr.rebx')
        sim.integrate(1.e4)
        sa = reboundx.Simulationarchive('binaries/gr.sa', 'binaries/gr.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)
    
    def test_gr_potential(self):
        sim = rebound.Simulation('binaries/twoplanets.bin')
        sim.save_to_file('binaries/gr_potential.sa', interval=1e3, delete_file=True)
        rebx = reboundx.Extras(sim)
        force = rebx.load_force('gr_potential')
        rebx.add_force(force)
        force.params['c'] = 1.e4
        ps = sim.particles
        rebx.save('binaries/gr_potential.rebx')
        sim.integrate(1.e4)
        sa = reboundx.Simulationarchive('binaries/gr_potential.sa', 'binaries/gr_potential.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)

    def test_radiation_forces(self):
        sim = rebound.Simulation('binaries/twoplanets.bin')
        sim.save_to_file('binaries/radiation_forces.sa', interval=1e3, delete_file=True)
        rebx = reboundx.Extras(sim)
        force = rebx.load_force('radiation_forces')
        rebx.add_force(force)
        force.params['c'] = 1.e4
        ps = sim.particles
        ps[1].params['beta'] = 0.5
        ps[2].params['beta'] = 0.3
        rebx.save('binaries/radiation_forces.rebx')
        sim.integrate(1.e4)
        sa = reboundx.Simulationarchive('binaries/radiation_forces.sa', 'binaries/radiation_forces.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)

    def test_central_force(self):
        sim = rebound.Simulation('binaries/twoplanets.bin')
        sim.save_to_file('binaries/central_force.sa', interval=1e3, delete_file=True)
        rebx = reboundx.Extras(sim)
        force = rebx.load_force('central_force')
        rebx.add_force(force)
        ps = sim.particles
        ps[0].params['Acentral'] = 1.e-3
        ps[0].params['gammacentral'] = -1
        rebx.save('binaries/central_force.rebx')
        sim.integrate(1.e4)
        sa = reboundx.Simulationarchive('binaries/central_force.sa', 'binaries/central_force.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)
    
    def test_gravitational_harmonics(self):
        sim = rebound.Simulation('binaries/twoplanets.bin')
        sim.save_to_file('binaries/gravitational_harmonics.sa', interval=1e3, delete_file=True)
        rebx = reboundx.Extras(sim)
        force = rebx.load_force('gravitational_harmonics')
        rebx.add_force(force)
        ps = sim.particles
        ps[0].params['J2'] = 1.e-3
        ps[0].params['J4'] = 1.e-3
        ps[0].params['R_eq'] = 1.e-3
        rebx.save('binaries/gravitational_harmonics.rebx')
        sim.integrate(1.e4)
        sa = reboundx.Simulationarchive('binaries/gravitational_harmonics.sa', 'binaries/gravitational_harmonics.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)

    def test_modify_mass(self):
        sim = rebound.Simulation('binaries/twoplanets.bin')
        sim.save_to_file('binaries/modify_mass.sa', interval=1e3, delete_file=True)
        rebx = reboundx.Extras(sim)
        mod = rebx.load_operator('modify_mass')
        rebx.add_operator(mod)
        ps = sim.particles
        ps[0].params['tau_mass'] = -1e4
        rebx.save('binaries/modify_mass.rebx')
        sim.integrate(1.e4)
        sa = reboundx.Simulationarchive('binaries/modify_mass.sa', 'binaries/modify_mass.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)


if __name__ == '__main__':
    unittest.main()
