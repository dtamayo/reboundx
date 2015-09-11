from ctypes import *
import os
import rebound
from rebound import Simulation

c_default = 10064.915
#Find the reboundx C library
pymodulespath = os.path.dirname(__file__)
try:
    clibreboundx = CDLL(pymodulespath + '/../libreboundx.so', RTLD_GLOBAL)
except:
    print("Cannot find library 'libreboundx.so'.")
    raise

class rebx_params_modify_orbits_forces(Structure):
    _fields_ = [("allocatedN", c_int),
                ("tau_a", POINTER(c_double)),
                ("tau_e", POINTER(c_double)),
                ("tau_inc", POINTER(c_double)),
                ("tau_omega", POINTER(c_double)),
                ("e_damping_p", c_double)]                     

class rebx_params_modify_orbits_direct(Structure):
    _fields_ = [("allocatedN", c_int),
                ("tau_a", POINTER(c_double)),
                ("tau_e", POINTER(c_double)),
                ("tau_inc", POINTER(c_double)),
                ("tau_omega", POINTER(c_double)),
                ("e_damping_p", c_double)]                     

class rebx_params_gr(Structure):
    _fields_ = [("c", c_double)]

class Extras(Structure):
    def __init__(self, sim):
        clibreboundx.rebx_initialize(byref(sim), byref(self))

    def __del__(self):
        if self._b_needsfree_ == 1:
            clibreboundx.rebx_free_pointers(byref(self))

    def add_modify_orbits_direct(self):
        clibreboundx.rebx_add_modify_orbits_direct(self.sim)
    
    def add_modify_orbits_forces(self):
        clibreboundx.rebx_add_modify_orbits_forces(self.sim)

    def check_c(self, c):
        if c is not None: # user passed c explicitly
            return c
      
        # c was not passed by user
         
        if self.sim.contents.G == 1: # if G = 1 (default) return default c
            return c_default
            
        u = self.sim.contents.units
        if not None in u.values(): # units are set
            from rebound import units
            c = units.convert_vel(c_default, 'au', 'yr2pi', u['length'], u['time'])
            return c
        else:
            raise ValueError("If you change G, you must pass c (speed of light) in appropriate units to add_gr, add_gr_potential, and add_gr_implicit.  Alternatively, set the units for the simulation.  See ipython_examples/GeneralRelativity.ipynb")
              
        return c_default
    def add_gr(self, c=None):
        c = self.check_c(c)
        clibreboundx.rebx_add_gr(self.sim, c_double(c))
    
    def add_gr_single_mass(self, c=None):
        c = self.check_c(c)
        clibreboundx.rebx_add_gr_single_mass(self.sim, c_double(c))

    def add_gr_potential(self, c=None):
        c = self.check_c(c)
        clibreboundx.rebx_add_gr_potential(self.sim, c_double(c))

'''
    @property
    def tau_a(self):
        return [self.params.contents.tau_a[i] for i in range(self.sim.contents.N)]
    @tau_a.setter
    def tau_a(self, tau_a):
        if len(tau_a) is not self.sim.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_a))(*tau_a)
        clibreboundx.rebx_set_tau_a(self.sim, byref(arr))

    @property
    def tau_e(self):
        return [self.params.contents.tau_e[i] for i in range(self.sim.contents.N)]
    @tau_e.setter
    def tau_e(self, tau_e):
        if len(tau_e) is not self.sim.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_e))(*tau_e)
        clibreboundx.rebx_set_tau_e(self.sim, byref(arr))

    @property
    def tau_inc(self):
        return [self.params.contents.tau_inc[i] for i in range(self.sim.contents.N)]
    @tau_inc.setter
    def tau_inc(self, tau_inc):
        if len(tau_inc) is not self.sim.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_inc))(*tau_inc)
        clibreboundx.rebx_set_tau_inc(self.sim, byref(arr))

    @property
    def tau_pomega(self):
        return [self.params.contents.tau_pomega[i] for i in range(self.sim.contents.N)]
    @tau_pomega.setter
    def tau_pomega(self, tau_pomega):
        if len(tau_pomega) is not self.sim.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_pomega))(*tau_pomega)
        clibreboundx.rebx_set_tau_pomega(self.sim, byref(arr))
'''
# Need to put fields after class definition because of self-referencing
Extras._fields_ = [("sim", POINTER(Simulation)),
                ("forces", POINTER(CFUNCTYPE(None, POINTER(Simulation)))),
                ("ptm", POINTER(CFUNCTYPE(None, POINTER(Simulation)))),
                ("Nforces", c_int),
                ("Nptm", c_int),
                ("modify_orbits_forces", rebx_params_modify_orbits_forces),
                ("modify_orbits_direct", rebx_params_modify_orbits_direct),
                ("gr", rebx_params_gr)]

'''
def set_e_damping_p(value):
    clibreboundx.rebx_set_e_damping_p(c_double(value))

def set_migration(tau_a):
    arr = (c_double * len(tau_a))(*tau_a)
    clibreboundx.rebx_set_migration(byref(arr),c_int(len(tau_a)))

def set_e_damping(tau_e, p=0.):
    c_double.in_dll(libreboundx, "e_damping_p").value = p
    arr = (c_double * len(tau_e))(*tau_e)
    clibreboundx.rebx_set_e_damping(byref(arr),c_int(len(tau_e)))
    
def set_e_damping(tau_e):
    arr = (c_double * len(tau_e))(*tau_e)
    clibreboundx.rebx_set_e_damping(byref(arr),c_int(len(tau_e)))
    
def set_i_damping(tau_i):
    arr = (c_double * len(tau_i))(*tau_i)
    clibreboundx.rebx_set_i_damping(byref(arr),c_int(len(tau_i)))

def set_peri_precession(tau_omega):
    arr = (c_double * len(tau_omega))(*tau_omega)
    clibreboundx.rebx_set_peri_precession(byref(arr), c_int(len(tau_omega)))
   
def set_peri_precession(gamma, Rc, podot):
    clibreboundx.set_peri_precession(c_double(gamma), c_double(Rc), c_double(podot))
'''
