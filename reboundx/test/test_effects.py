import rebound
import reboundx
import unittest

class TestForces(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(a=1., e=0.2)
        self.rebx = reboundx.Extras(self.sim)

    def test_loadforce(self):
        gr = self.rebx.load_force('gr')
    
    def test_loadforcewrongname(self):
        with self.assertRaises(RuntimeError):
            gr = self.rebx.load_force('gr2')

    def test_customforce(self):
        cust = self.rebx.create_force('myforce')
        def myforce(sim, force, particles, N):
            sim.contents.particles[1].ax += 1.e-4
        cust.update_accelerations = myforce
        cust.force_type = 'pos'
        self.rebx.add_force(cust)
        self.sim.integrate(10)
        self.assertGreater(self.sim.particles[1].pomega, -0.01)

    def test_customnoforce(self):
        cust = self.rebx.create_force('myforce')
        cust.force_type = 'pos'
        with self.assertRaises(RuntimeError):
            self.rebx.add_force(cust)

    def test_customnoforcetype(self):
        cust = self.rebx.create_force('myforce')
        def myforce(sim, force, particles, N):
            sim.contents.particles[1].ax += 1.e-4
        cust.update_accelerations = myforce
        with self.assertRaises(RuntimeError):
            self.rebx.add_force(cust)
    
    def test_getforce(self):
        gr = self.rebx.load_force('gr')
        gr.params['gr_source'] = 3
        newgr = self.rebx.get_force('gr')
        self.assertEqual(newgr.params['gr_source'], 3)
    
    def test_getforcenotfirst(self):
        gr = self.rebx.load_force('gr')
        gr.params['gr_source'] = 3
        cust = self.rebx.create_force('myforce')
        newgr = self.rebx.get_force('gr')
        self.assertEqual(newgr.params['gr_source'], 3)

    def test_forcenotfound(self):
        with self.assertRaises(AttributeError):
            gr = self.rebx.get_force('gr')
    
    def test_removeforce(self):
        gr = self.rebx.load_force('gr')
        self.rebx.add_force(gr)
        self.rebx.remove_force(gr)
    
    def test_removeforcenotfirst(self):
        gr = self.rebx.load_force('gr')
        self.rebx.add_force(gr)
        cust = self.rebx.create_force('myforce')
        def myforce(sim, force, particles, N):
            sim.contents.particles[1].ax += 1.e-4
        cust.update_accelerations = myforce
        cust.force_type = 'pos'
        self.rebx.add_force(cust)
        self.rebx.remove_force(gr)
    
    def test_removenonforce(self):
        with self.assertRaises(TypeError):
            self.rebx.remove_force(self.sim)

    def test_removeforcenotattached(self):
        gr = self.rebx.load_force('gr')
        with self.assertRaises(AttributeError):
            self.rebx.remove_force(gr)

class TestOperators(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(a=1., e=0.2)
        self.rebx = reboundx.Extras(self.sim)

    def test_loadoperator(self):
        mm = self.rebx.load_operator('modify_mass')
    
    def test_loadoperatorwrongname(self):
        with self.assertRaises(RuntimeError):
            mm = self.rebx.load_force('modify_mass2')

    def test_customoperator(self):
        cust = self.rebx.create_operator('myoperator')
        def mystep(sim, operator, dt):
            sim.contents.particles[1].x += 1.e-4
        cust.step_function = mystep
        cust.operator_type = 'updater'
        self.rebx.add_operator(cust)
        self.sim.integrate(10)
        self.assertGreater(self.sim.particles[1].pomega, 1.e-3)

    def test_customopnostep(self):
        cust = self.rebx.create_operator('myoperator')
        cust.operator_type = 'updater'
        with self.assertRaises(RuntimeError):
            self.rebx.add_operator(cust)

    def test_customopnooperatortype(self):
        cust = self.rebx.create_operator('myoperator')
        def mystep(sim, operator, dt):
            sim.contents.particles[1].x += 1.e-4
        cust.step_function = mystep
        with self.assertRaises(RuntimeError):
            self.rebx.add_operator(cust)

    def test_getoperator(self):
        mm = self.rebx.load_operator('modify_mass')
        mm.params['gr_source'] = 3
        newmm = self.rebx.get_operator('modify_mass')
        self.assertEqual(newmm.params['gr_source'], 3)
    
    def test_getoperatornotfirst(self):
        mm = self.rebx.load_operator('modify_mass')
        cust = self.rebx.create_operator('myoperator')
        mm.params['gr_source'] = 3
        newmm = self.rebx.get_operator('modify_mass')
        self.assertEqual(newmm.params['gr_source'], 3)

    def test_operatornotfound(self):
        with self.assertRaises(AttributeError):
            mm = self.rebx.get_operator('modify_mass')
    
    def test_removeoperator(self):
        mm = self.rebx.load_operator('modify_mass')
        self.rebx.add_operator(mm)
        self.rebx.remove_operator(mm)

    def test_removeoperatorprepost(self):
        self.sim.integrator = "whfast"
        mm = self.rebx.load_operator('modify_mass')
        self.rebx.add_operator(mm)
        self.rebx.remove_operator(mm)
    
    def test_removemanyoperatorsteps(self):
        mm = self.rebx.load_operator('modify_mass')
        cust = self.rebx.create_operator('myoperator')
        def mystep(sim, operator, dt):
            sim.contents.particles[1].x += 1.e-4
        cust.step_function = mystep
        cust.operator_type = 'updater'
        self.rebx.add_operator(mm, dtfraction=0.5, timing='post')
        self.rebx.add_operator(cust, dtfraction=0.5, timing='post')
        self.rebx.add_operator(mm, dtfraction=0.5, timing='post')
        self.rebx.add_operator(cust, dtfraction=0.5, timing='post')
        self.rebx.add_operator(mm, dtfraction=0.5, timing='post')
        self.rebx.add_operator(cust, dtfraction=0.5, timing='pre')
        self.rebx.add_operator(mm, dtfraction=0.5, timing='pre')
        self.rebx.add_operator(cust, dtfraction=0.5, timing='pre')
        self.rebx.add_operator(mm, dtfraction=0.5, timing='pre')
        self.rebx.add_operator(cust, dtfraction=0.5, timing='pre')
        self.rebx.remove_operator(mm)
    
    def test_removenonoperator(self):
        with self.assertRaises(TypeError):
            self.rebx.remove_operator(self.sim)

    def test_removeoperatornotattached(self):
        mm = self.rebx.load_operator('modify_mass')
        with self.assertRaises(AttributeError):
            self.rebx.remove_operator(mm)

class TestAddOperator(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(a=1., e=0.2)
        self.rebx = reboundx.Extras(self.sim)
        self.cust = self.rebx.create_operator('addop')
        self.rebx.register_param('ctr', 'REBX_TYPE_INT')
        self.cust.params['ctr'] = 0
        def mystep(sim, operator, dt):
            operator.contents.params['ctr'] += 1
        self.cust.step_function = mystep
    
    def test_ias15updater(self):
        self.sim.integrator='ias15'
        self.cust.operator_type = 'updater'
        self.rebx.add_operator(self.cust)
        self.sim.step()
        self.assertEqual(self.cust.params['ctr'], 1)
    
    def test_ias15recorder(self):
        self.sim.integrator='ias15'
        self.cust.operator_type = 'recorder'
        self.rebx.add_operator(self.cust)
        self.sim.step()
        self.assertEqual(self.cust.params['ctr'], 1)
    
    def test_whfastupdater(self):
        self.sim.integrator='whfast'
        self.cust.operator_type = 'updater'
        self.rebx.add_operator(self.cust)
        self.sim.step()
        self.assertEqual(self.cust.params['ctr'], 2) # whfast = sec order scheme
    
    def test_whfastrecorder(self):
        self.sim.integrator='whfast'
        self.cust.operator_type = 'recorder'
        self.rebx.add_operator(self.cust)
        self.sim.step()
        self.assertEqual(self.cust.params['ctr'], 1)

    def test_mercuriusupdater(self):
        self.sim.integrator='mercurius'
        self.cust.operator_type = 'updater'
        with self.assertRaises(RuntimeError):
            self.rebx.add_operator(self.cust)
    
    def test_mercuriusrecorder(self):
        self.sim.integrator='mercurius'
        self.cust.operator_type = 'recorder'
        self.rebx.add_operator(self.cust)
        self.sim.step()
        self.assertEqual(self.cust.params['ctr'], 1)

if __name__ == '__main__':
    unittest.main()

