import rebound
#import reboundx
#from reboundx import data
import unittest
import math
import numpy as np
from ctypes import c_uint, c_uint8, c_uint32, c_uint64

def mycomp(obj1, obj2):
    if type(obj1) != type(obj2):
        return False
    for attr in [attr for attr in dir(obj1) if not attr.startswith('_')]:
        if getattr(obj1, attr) != getattr(obj2, attr):
            return False
    return True

class TestEffectParams(unittest.TestCase):
    def setUp(self):
        pass#self.sim = rebound.Simulation()
        '''
        data.add_earths(self.sim, ei=1.e-3) 
        self.rebx = reboundx.Extras(self.sim)
        self.gr = self.rebx.add("gr")
        self.gr.params["a"] = 1.2
        self.gr.params["b"] = 1.7
        self.gr.params["N"] = 14
        '''

    def tearDown(self):
        self.sim = None

'''
    def test_access(self):
        self.assertAlmostEqual(self.gr.params["a"], 1.2, delta=1.e-15)
        self.assertAlmostEqual(self.gr.params["b"], 1.7, delta=1.e-15)
        self.assertEqual(self.gr.params["N"], 14)

    def test_types(self):
        for t in [int, float, np.int32, np.float64]:
            var = t(3) # make instance of type
            self.gr.params[t.__name__] = var
            self.assertAlmostEqual(self.gr.params[t.__name__], var)
            self.assertEqual(type(self.gr.params[t.__name__]), type(var))
        for t in [c_uint, c_uint32]:
            var = t(3) # make instance of type
            self.gr.params[t.__name__] = var
            self.assertAlmostEqual(self.gr.params[t.__name__].value, var.value)
            self.assertEqual(type(self.gr.params[t.__name__]), type(var))
        self.gr.params["orbit"] = rebound.Orbit()
        self.assertEqual(type(self.gr.params["orbit"]), rebound.Orbit)

    def test_iter(self):
        with self.assertRaises(AttributeError):
            for p in self.gr.params:
                pass

    def test_length(self):
        self.assertEqual(len(self.gr.params), 3)

    def test_del(self):
        del self.gr.params["b"]
        with self.assertRaises(AttributeError):
            self.gr.params["b"]
        self.assertEqual(len(self.gr.params), 2)

        del self.gr.params["a"]
        with self.assertRaises(AttributeError):
            self.gr.params["a"]
        self.assertEqual(len(self.gr.params), 1)

        del self.gr.params["N"]
        with self.assertRaises(AttributeError):
            self.gr.params["N"]
        self.assertEqual(len(self.gr.params), 0)
'''

class TestParticleParams(unittest.TestCase):
    def setUp(self):
        pass#self.sim = rebound.Simulation()
        '''
        self.rebx = reboundx.Extras(self.sim)
        data.add_earths(self.sim, ei=1.e-3) 
        self.sim.particles[0].params["a"] = 1.2
        self.sim.particles[0].params["b"] = 1.7
        self.sim.particles[0].params["N"] = 14
        '''
    def tearDown(self):
        self.sim = None
'''    
    def test_access(self):
        self.assertAlmostEqual(self.sim.particles[0].params["a"], 1.2, delta=1.e-15)
        self.assertAlmostEqual(self.sim.particles[0].params["b"], 1.7, delta=1.e-15)
        self.assertEqual(self.sim.particles[0].params["N"], 14)

    def test_iter(self):
        with self.assertRaises(AttributeError):
            for p in self.sim.particles[0].params:
                pass

    def test_update(self):
        self.sim.particles[0].params["a"] = 42.
        self.assertAlmostEqual(self.sim.particles[0].params["a"], 42., delta=1.e-15)

    def test_update_with_wrong_type(self):
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["N"] = 37.2
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["a"] = 37
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["a"] = [1., 2., 3.]
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["a"] = np.array([1., 2., 3.])


    def test_length(self):
        self.assertEqual(len(self.sim.particles[0].params), 3)

    def test_del(self):
        del self.sim.particles[0].params["b"]
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["b"]
        self.assertEqual(len(self.sim.particles[0].params), 2)

        del self.sim.particles[0].params["a"]
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["a"]
        self.assertEqual(len(self.sim.particles[0].params), 1)

        del self.sim.particles[0].params["N"]
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["N"]
        self.assertEqual(len(self.sim.particles[0].params), 0)
'''

class TestRebxNotAttached(unittest.TestCase):
    def test_rebx_not_attached(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        #with self.assertRaises(AttributeError):
        #    self.sim.particles[0].params["a"] = 7

if __name__ == '__main__':
    unittest.main()
