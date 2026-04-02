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
        
from ctypes import Structure, c_double, cast, POINTER
class Mystruct(Structure):
    _fields_ = [('dt', c_double),
                ('c', c_double)]

class TestParams(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(a=1.)
        self.rebx = reboundx.Extras(self.sim)
        self.gr = self.rebx.load_force("gr")
        self.mm = self.rebx.load_operator("modify_mass")
        self.p = self.sim.particles[1]

    def tearDown(self):
        self.sim = None

    def test_adddouble(self):
        self.gr.params['c'] = 1.2
        self.mm.params['c'] = 1.4
        self.p.params['c'] = 1.7
        self.assertAlmostEqual(self.gr.params["c"], 1.2, delta=1.e-15)
        self.assertAlmostEqual(self.mm.params["c"], 1.4, delta=1.e-15)
        self.assertAlmostEqual(self.p.params["c"], 1.7, delta=1.e-15)

    def test_updatedouble(self):
        self.gr.params['c'] = 1.2
        self.mm.params['c'] = 1.4
        self.p.params['c'] = 1.7
        
        self.gr.params['c'] = 2.2
        self.mm.params['c'] = 2.4
        self.p.params['c'] = 2.7
        self.assertAlmostEqual(self.gr.params["c"], 2.2, delta=1.e-15)
        self.assertAlmostEqual(self.mm.params["c"], 2.4, delta=1.e-15)
        self.assertAlmostEqual(self.p.params["c"], 2.7, delta=1.e-15)
    
    def test_updatedoublecopy(self):
        a = 1.2
        self.gr.params['c'] = a

        a = 3.7
        self.assertAlmostEqual(self.gr.params["c"], 1.2, delta=1.e-15) # shouldn't reflect update because copied
    
    def test_addint(self):
        self.gr.params['gr_source'] = 1
        self.mm.params['gr_source'] = 2
        self.p.params['gr_source'] = 3
        self.assertAlmostEqual(self.gr.params["gr_source"], 1, delta=1.e-15)
        self.assertAlmostEqual(self.mm.params["gr_source"], 2, delta=1.e-15)
        self.assertAlmostEqual(self.p.params["gr_source"], 3, delta=1.e-15)
    
    def test_updateint(self):
        self.gr.params['gr_source'] = 1
        self.mm.params['gr_source'] = 2 
        self.p.params['gr_source'] = 3 
        
        self.gr.params['gr_source'] = 11 
        self.mm.params['gr_source'] = 22
        self.p.params['gr_source'] = 33
        self.assertEqual(self.gr.params["gr_source"], 11)
        self.assertEqual(self.mm.params["gr_source"], 22)
        self.assertEqual(self.p.params["gr_source"], 33)
   
    def test_updateintcopy(self):
        a = 1
        self.gr.params['gr_source'] = a

        a = 3
        self.assertEqual(self.gr.params["gr_source"], 1) # shouldn't reflect update because copied
    
    def test_addforce(self):
        self.gr.params['c'] = 3.5
        self.p.params['force'] = self.gr
        newgr = self.p.params['force']
        self.assertAlmostEqual(newgr.params['c'], 3.5, delta=1.e-15)

    def test_addnonforce(self):
        with self.assertRaises(AttributeError):
            self.p.params['force'] = 3.

    def test_updateforceptr(self):
        self.gr.params['c'] = 3.5
        self.p.params['force'] = self.gr
        self.gr.params['c'] = 4.5
        newgr = self.p.params['force']
        self.assertAlmostEqual(newgr.params['c'], 4.5, delta=1.e-15)
    
    def test_updateforce(self):
        self.p.params['force'] = self.gr
        self.gr.params['c'] = 3.5

        newforce = self.rebx.load_force('gr')
        newforce.params['c'] = -3.5
        self.p.params['force'] = newforce
        newgr = self.p.params['force']
        self.assertAlmostEqual(newgr.params['c'], -3.5, delta=1.e-15)
    
    def test_addcustomstruct(self):
        self.rebx.register_param('my_new_struct', 'REBX_TYPE_POINTER')
        s = Mystruct()
        s.dt = 0.1
        s.c = 3.5
        self.gr.params['my_new_struct'] = s
        new_s = self.gr.params['my_new_struct']
        new_s = cast(new_s, POINTER(Mystruct)).contents
        self.assertAlmostEqual(new_s.dt, 0.1, delta=1.e-15)
        self.assertAlmostEqual(new_s.c, 3.5, delta=1.e-15)

    def test_updatecustomstructptr(self):
        self.rebx.register_param('my_new_struct', 'REBX_TYPE_POINTER')
        s = Mystruct()
        s.dt = 0.1
        s.c = 3.5
        self.gr.params['my_new_struct'] = s

        s.dt = 1.1
        new_s = self.gr.params['my_new_struct']
        new_s = cast(new_s, POINTER(Mystruct)).contents
        self.assertAlmostEqual(new_s.dt, 1.1, delta=1.e-15)
    
    def test_updatecustomstruct(self):
        self.rebx.register_param('my_new_struct', 'REBX_TYPE_POINTER')
        s = Mystruct()
        s.dt = 0.1
        s.c = 3.5
        self.gr.params['my_new_struct'] = s

        s2 = Mystruct()
        s2.dt = -0.1
        s2.c = -3.5
        self.gr.params['my_new_struct'] = s2

        new_s = self.gr.params['my_new_struct']
        new_s = cast(new_s, POINTER(Mystruct)).contents
        self.assertAlmostEqual(new_s.dt, -0.1, delta=1.e-15)
        self.assertAlmostEqual(new_s.c, -3.5, delta=1.e-15)

    def test_getnotregistered(self):
        with self.assertRaises(AttributeError):
            b = self.gr.params['asd;flkj']

    def test_notregistered(self):
        with self.assertRaises(AttributeError):
            self.gr.params['asldkjf'] = 1.2
    
    def test_registerexisting(self):
        with self.assertRaises(RuntimeError):
            self.rebx.register_param('c', 'REBX_TYPE_INT')

    def test_not_attached(self):
        with self.assertRaises(AttributeError):
            b = self.gr.params['beta']
    
    def test_custom_not_attached(self):
        self.rebx.register_param('my_custom', 'REBX_TYPE_POINTER')
        with self.assertRaises(AttributeError):
            b = self.gr.params['my_custom']

    def test_newdouble(self):
        self.rebx.register_param('my_new_double', 'REBX_TYPE_DOUBLE')
        self.gr.params['my_new_double'] = 1.2
        self.assertAlmostEqual(self.gr.params["my_new_double"], 1.2, delta=1.e-15)
    
    def test_newint(self):
        self.rebx.register_param('my_new_int', 'REBX_TYPE_INT')
        self.gr.params['my_new_int'] = 2
        self.assertEqual(self.gr.params["my_new_int"], 2)

    def test_length(self):
        self.gr.params['c'] = 1.3
        self.gr.params['gr_source'] = 7
        self.gr.params['tau_mass'] = 3.2
        self.assertEqual(len(self.gr.params), 3)

    def test_iter(self):
        with self.assertRaises(AttributeError):
            for p in self.gr.params:
                pass

    def test_del(self):
        with self.assertRaises(AttributeError):
            del self.gr.params["b"]

if __name__ == '__main__':
    unittest.main()
