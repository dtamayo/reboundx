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

class TestSegFaults(unittest.TestCase):
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
