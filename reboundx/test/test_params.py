import rebound
import reboundx
from reboundx import data
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
        self.sim = rebound.Simulation()
        data.add_earths(self.sim, ei=1.e-3) 
        self.rebx = reboundx.Extras(self.sim)
        self.gr = self.rebx.add("gr")
        self.gr.params["a"] = 1.2
        self.gr.params["b"] = 1.7
        self.gr.params["N"] = 14

    def tearDown(self):
        self.sim = None
    
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

class TestParticleParams(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.rebx = reboundx.Extras(self.sim)
        data.add_earths(self.sim, ei=1.e-3) 
        self.sim.particles[0].params["a"] = 1.2
        self.sim.particles[0].params["b"] = 1.7
        self.sim.particles[0].params["N"] = 14

    def tearDown(self):
        self.sim = None
    
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

class TestArrays(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.rebx = reboundx.Extras(self.sim)
        data.add_earths(self.sim, ei=1.e-3) 
        self.array=np.array([1.,2.,3.,4.])
        self.intarray=np.array([1,2,3,4])
        self.cuintarray=np.array([c_uint(1),c_uint(2),c_uint(3),c_uint(4)], dtype=object)
        self.cuint32array=np.array([c_uint32(1),c_uint32(2),c_uint32(3),c_uint32(4)], dtype=object)
        self.ndarray = np.array([[[1.,2.,3.],[4.,5.,6.]], [[1.,2.,3.],[4.,5.,6.]],[[1.,2.,3.],[4.,5.,6.]], [[1.,2.,3.],[4.,5.,6.]]])
        self.orbitarray =np.array(self.sim.calculate_orbits(), dtype=object) 
        o1 = self.sim.particles[1].orbit
        o2 = self.sim.particles[2].orbit
        self.objndarray = np.array([[[o1,o2],[o2,o1],[o1,o1]], [[o1,o2],[o2,o1],[o1,o1]],[[o1,o2],[o2,o1],[o1,o1]],[[o1,o2],[o2,o1],[o1,o1]]], dtype=object)
        
        self.sim.particles[0].params["array"] = self.array
        self.sim.particles[0].params["intarray"] = self.intarray
        self.sim.particles[0].params["cuintarray"] = self.cuintarray
        self.sim.particles[0].params["cuint32array"] = self.cuint32array 
        self.sim.particles[1].params["ndarray"] = self.ndarray
        self.sim.particles[0].params["orbitarray"] = self.orbitarray
        self.sim.particles[1].params["objndarray"] = self.objndarray
        self.sim.particles[1].params["scalar"] = 3

    def tearDown(self):
        self.sim = None
   
    def test_objectarray(self): # param is array of custom classes with dtype=object
        for i, orbit in enumerate(self.sim.particles[0].params["orbitarray"]):
            self.assertEqual(mycomp(orbit, self.sim.particles[i+1].orbit), True)
        
    '''
    Seems like with some versions of numpy, making an array of objects defaults to setting dtype=object (one on travis), and others makes flexible dtype (mine), so don't test.
    def test_incorrectarray(self): # param is array of custom classes but dtype=object not set
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["badarray"] = np.array(self.sim.calculate_orbits())
    '''

    def test_arrays(self): # 3d array of floats
        self.assertEqual(np.array_equal(self.array, self.sim.particles[0].params["array"]), True)
        self.assertEqual(np.array_equal(self.intarray, self.sim.particles[0].params["intarray"]), True)
        self.assertEqual(np.array_equal(self.ndarray, self.sim.particles[1].params["ndarray"]), True)
        self.assertAlmostEqual(self.sim.particles[0].params["cuintarray"][0].value, self.cuintarray[0].value)
        self.assertAlmostEqual(self.sim.particles[0].params["cuint32array"][0].value, self.cuint32array[0].value)

    def test_objndarray(self): # 3d array of rebound.Orbits
        flat = self.objndarray.flatten()
        for i, orbit in enumerate(self.sim.particles[1].params["objndarray"].flatten()):
            self.assertEqual(mycomp(orbit, flat[i]), True)
    
    def test_shared_memory_objects(self):
        q = self.sim.particles[0].params["orbitarray"]
        q[0].a = 298
        self.assertAlmostEqual(self.sim.particles[0].params["orbitarray"][0].a, 298, delta=1.e-15)

    def test_diff_shape(self):
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["array"] = np.array([1.,2.])
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["list"] = [1.,2.]
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["orbitarray"] = np.array([self.sim.particles[1].orbit])
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["orbitlist"] = [self.sim.particles[1].orbit]
        with self.assertRaises(AttributeError):
            self.sim.particles[1].params["scalar"] = [1.,2.]
        with self.assertRaises(AttributeError):
            self.sim.particles[1].params["scalar"] = np.array([1.,2.])
        with self.assertRaises(AttributeError):
            self.sim.particles[1].params["scalar"] = [1,2]
        with self.assertRaises(AttributeError):
            self.sim.particles[1].params["scalar"] = np.array([1,2])

    def test_update_array(self):
        newlist = [4.,3.,2.,1.]
        newarray = np.array(newlist)
        self.sim.particles[0].params["array"] = newarray
        self.assertEqual(np.array_equal(self.sim.particles[0].params["array"], newarray), True)
    
    def test_iter(self):
        with self.assertRaises(AttributeError):
            for p in self.sim.particles[0].params:
                pass

    def test_length(self):
        self.assertEqual(len(self.sim.particles[0].params), 5)

    def test_del(self):
        del self.sim.particles[0].params["array"]
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["array"]
        self.assertEqual(len(self.sim.particles[0].params), 4)

        del self.sim.particles[0].params["orbitarray"]
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["orbitarray"]
        self.assertEqual(len(self.sim.particles[0].params), 3)

class TestRebxNotAttached(unittest.TestCase):
    def test_rebx_not_attached(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        with self.assertRaises(AttributeError):
            self.sim.particles[0].params["a"] = 7

if __name__ == '__main__':
    unittest.main()
