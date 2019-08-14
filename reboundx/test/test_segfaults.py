import rebound
import reboundx
import unittest

def make():
    sim = rebound.Simulation()
    sim.add(m=1.)
    sim.add(a=0.1, e=0.1)
    rebx = reboundx.Extras(sim)
    return sim, rebx
        
def makegr():
    sim, rebx = make()
    gr = rebx.load_force('gr')
    rebx.add_force(gr)
    gr.params['c'] = 1e2
    return sim, rebx

def dropsim():
    sim, rebx = makegr()
    return rebx

def droprebx():
    sim, rebx = makegr()
    return sim

def testEffect(reb_sim):
    reb_sim.contents.particles[0].params['c'] = 1.
    
class TestSegFaults(unittest.TestCase):
    def test_detach_preserves_force(self):
        sim, rebx = make()
        sim.additional_forces = testEffect
        addmass = rebx.load_operator('modify_mass')
        rebx.add_operator(addmass, dtfraction=1.,  timing="pre")
        rebx.detach(sim)
        sim.step()
        self.assertEqual(sim.particles[0].params['c'], 1)
    
    def test_detach_preserves_pretm(self):
        sim, rebx = make()
        sim.pre_timestep_modifications = testEffect
        addmass = rebx.load_operator('modify_mass')
        rebx.add_operator(addmass, dtfraction=1.,  timing="post")
        rebx.detach(sim)
        sim.step()
        self.assertEqual(sim.particles[0].params['c'], 1)
    
    def test_detach_preserves_pretm(self):
        sim, rebx = make()
        sim.post_timestep_modifications = testEffect
        gr = rebx.load_force('gr')
        rebx.add_force(gr)
        rebx.detach(sim)
        sim.step()
        self.assertEqual(sim.particles[0].params['c'], 1)

    def test_detach(self):
        sim, rebx = makegr()
        rebx.detach(sim)
        sim.integrate(10)
        self.assertLess(sim.particles[1].pomega, 1.e-10) # should be 0

    def test_create_in_function(self):
        sim, rebx = makegr()
        sim.integrate(10)
        self.assertGreater(sim.particles[1].pomega, 0.01)

    def test_droprebx(self):
        sim = droprebx()
        sim.integrate(10)
        self.assertGreater(sim.particles[1].pomega, 0.01)
    
    def test_dropsim(self):
        rebx = dropsim()
        with self.assertRaises(AttributeError):
            rebx.load_force("gr")

    def test_delete(self):
        for i in range(10):
            sim, rebx = makegr()
            sim.integrate(10)
        self.assertGreater(sim.particles[1].pomega, 0.01)

    def test_rebx_not_attached(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        with self.assertRaises(AttributeError):
            sim.particles[0].params["a"] = 7

if __name__ == '__main__':
    unittest.main()
