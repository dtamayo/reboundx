import rebound
import reboundx
import unittest
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class TestMachineIndependent(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(a=1., e=0.2)
        self.rebx = reboundx.Extras(self.sim)

    def test_modify_orbits_forces(self):
        name = "modify_orbits_forces"
        try:
            rebbin = os.path.join(THIS_DIR, 'binaries/'+name+'.sa')
            rebxbin = os.path.join(THIS_DIR, 'binaries/'+name+'.rebx')
            sa = reboundx.SimulationArchive(rebbin, rebxbin)
        except:
            try:
                sim = rebound.Simulation('binaries/twoplanets.bin')
            except:
                sim = rebound.Simulation()
                sim.add(m=1.)
                sim.add(m=1.e-4, a=1., e=0.1)
                sim.add(m=1.e-4, a=2, e=0.1, inc=0.2)
                sim.integrator="whfast"
                sim.dt = sim.particles[1].P/100
                sim.move_to_com()
                sim.save('binaries/twoplanets.bin')
           
            sim.automateSimulationArchive('binaries/modify_orbits_forces.sa', interval=1e3, deletefile=True)
            rebx = reboundx.Extras(sim)
            mod = rebx.load_force('binaries/modify_orbits_forces')
            rebx.add_force(mod)
            ps = sim.particles
            ps[1].params['tau_a'] = -1e4
            ps[1].params['tau_e'] = -1e3
            ps[2].params['tau_e'] = -1e3
            ps[2].params['tau_inc'] = -1e3
            rebx.save('binaries/modify_orbits_forces.rebx')
            sim.integrate(1.e4)
            sa = reboundx.SimulationArchive('binaries/modify_orbits_forces.sa', 'binaries/modify_orbits_forces.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)

    def test_gr_full(self):
        name = 'gr_full'
        try:
            rebbin = os.path.join(THIS_DIR, 'binaries/'+name+'.sa')
            rebxbin = os.path.join(THIS_DIR, 'binaries/'+name+'.rebx')
            sa = reboundx.SimulationArchive(rebbin, rebxbin)
        except:
            sim = rebound.Simulation('binaries/twoplanets.bin')
            sim.automateSimulationArchive('binaries/'+name+'.sa', interval=1e3, deletefile=True)
            rebx = reboundx.Extras(sim)
            force = rebx.load_force(name)
            rebx.add_force(force)
            force.params['c'] = 1.e4
            ps = sim.particles
            rebx.save('binaries/'+name+'.rebx')
            sim.integrate(1.e4)
            sa = reboundx.SimulationArchive('binaries/'+name+'.sa', 'binaries/'+name+'.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)
    
    def test_gr(self):
        name = 'gr'
        try:
            rebbin = os.path.join(THIS_DIR, 'binaries/'+name+'.sa')
            rebxbin = os.path.join(THIS_DIR, 'binaries/'+name+'.rebx')
            sa = reboundx.SimulationArchive(rebbin, rebxbin)
        except:
            sim = rebound.Simulation('binaries/twoplanets.bin')
            sim.automateSimulationArchive('binaries/'+name+'.sa', interval=1e3, deletefile=True)
            rebx = reboundx.Extras(sim)
            force = rebx.load_force(name)
            rebx.add_force(force)
            force.params['c'] = 1.e4
            ps = sim.particles
            rebx.save('binaries/'+name+'.rebx')
            sim.integrate(1.e4)
            sa = reboundx.SimulationArchive('binaries/'+name+'.sa', 'binaries/'+name+'.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)
    
    def test_gr_potential(self):
        name = 'gr_potential'
        try:
            rebbin = os.path.join(THIS_DIR, 'binaries/'+name+'.sa')
            rebxbin = os.path.join(THIS_DIR, 'binaries/'+name+'.rebx')
            sa = reboundx.SimulationArchive(rebbin, rebxbin)
        except:
            sim = rebound.Simulation('binaries/twoplanets.bin')
            sim.automateSimulationArchive('binaries/'+name+'.sa', interval=1e3, deletefile=True)
            rebx = reboundx.Extras(sim)
            force = rebx.load_force(name)
            rebx.add_force(force)
            force.params['c'] = 1.e4
            ps = sim.particles
            rebx.save('binaries/'+name+'.rebx')
            sim.integrate(1.e4)
            sa = reboundx.SimulationArchive('binaries/'+name+'.sa', 'binaries/'+name+'.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)

    def test_radiation_forces(self):
        name = 'radiation_forces'
        try:
            rebbin = os.path.join(THIS_DIR, 'binaries/'+name+'.sa')
            rebxbin = os.path.join(THIS_DIR, 'binaries/'+name+'.rebx')
            sa = reboundx.SimulationArchive(rebbin, rebxbin)
        except:
            sim = rebound.Simulation('binaries/twoplanets.bin')
            sim.automateSimulationArchive('binaries/'+name+'.sa', interval=1e3, deletefile=True)
            rebx = reboundx.Extras(sim)
            force = rebx.load_force(name)
            rebx.add_force(force)
            force.params['c'] = 1.e4
            ps = sim.particles
            ps[1].params['beta'] = 0.5
            ps[2].params['beta'] = 0.3
            rebx.save('binaries/'+name+'.rebx')
            sim.integrate(1.e4)
            sa = reboundx.SimulationArchive('binaries/'+name+'.sa', 'binaries/'+name+'.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)

    def test_central_force(self):
        name = 'central_force'
        try:
            rebbin = os.path.join(THIS_DIR, 'binaries/'+name+'.sa')
            rebxbin = os.path.join(THIS_DIR, 'binaries/'+name+'.rebx')
            sa = reboundx.SimulationArchive(rebbin, rebxbin)
        except:
            sim = rebound.Simulation('binaries/twoplanets.bin')
            sim.automateSimulationArchive('binaries/'+name+'.sa', interval=1e3, deletefile=True)
            rebx = reboundx.Extras(sim)
            force = rebx.load_force(name)
            rebx.add_force(force)
            ps = sim.particles
            ps[0].params['Acentral'] = 1.e-3
            ps[0].params['gammacentral'] = -1
            rebx.save('binaries/'+name+'.rebx')
            sim.integrate(1.e4)
            sa = reboundx.SimulationArchive('binaries/'+name+'.sa', 'binaries/'+name+'.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)
    
    def test_gravitational_harmonics(self):
        name = 'gravitational_harmonics'
        try:
            rebbin = os.path.join(THIS_DIR, 'binaries/'+name+'.sa')
            rebxbin = os.path.join(THIS_DIR, 'binaries/'+name+'.rebx')
            sa = reboundx.SimulationArchive(rebbin, rebxbin)
        except:
            sim = rebound.Simulation('binaries/twoplanets.bin')
            sim.automateSimulationArchive('binaries/'+name+'.sa', interval=1e3, deletefile=True)
            rebx = reboundx.Extras(sim)
            force = rebx.load_force(name)
            rebx.add_force(force)
            ps = sim.particles
            ps[0].params['J2'] = 1.e-3
            ps[0].params['J4'] = 1.e-3
            ps[0].params['R_eq'] = 1.e-3
            rebx.save('binaries/'+name+'.rebx')
            sim.integrate(1.e4)
            sa = reboundx.SimulationArchive('binaries/'+name+'.sa', 'binaries/'+name+'.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)

    def test_modify_mass(self):
        name = 'modify_mass'
        try:
            rebbin = os.path.join(THIS_DIR, 'binaries/'+name+'.sa')
            rebxbin = os.path.join(THIS_DIR, 'binaries/'+name+'.rebx')
            sa = reboundx.SimulationArchive(rebbin, rebxbin)
        except:
            sim = rebound.Simulation('binaries/twoplanets.bin')
            sim.automateSimulationArchive('binaries/'+name+'.sa', interval=1e3, deletefile=True)
            rebx = reboundx.Extras(sim)
            mod = rebx.load_operator(name)
            rebx.add_operator(mod)
            ps = sim.particles
            ps[0].params['tau_mass'] = -1e4
            rebx.save('binaries/'+name+'.rebx')
            sim.integrate(1.e4)
            sa = reboundx.SimulationArchive('binaries/'+name+'.sa', 'binaries/'+name+'.rebx')

        simf, rebx = sa[-1]
        sim,  rebx = sa[0]
        sim.integrate(simf.t)
        self.assertEqual(sim.particles[0].x, simf.particles[0].x)


if __name__ == '__main__':
    unittest.main()
