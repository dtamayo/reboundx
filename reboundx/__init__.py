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
        clibreboundx.rebx_initialize(byref(sim), byref(self)) # Use memory address ctypes allocated for rebx Structure in C
        self.add_Particle_props()

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

    def add_Particle_props(self):
        @property
        def tau_a(self):
            clibreboundx.rebx_get_tau_a.restype = c_double
            return clibreboundx.rebx_get_tau_a(byref(self))
        @tau_a.setter
        def tau_a(self, value):
            clibreboundx.rebx_set_tau_a(byref(self), c_double(value))
        @property
        def tau_e(self):
            clibreboundx.rebx_get_tau_e.restype = c_double
            return clibreboundx.rebx_get_tau_e(byref(self))
        @tau_e.setter
        def tau_e(self, value):
            clibreboundx.rebx_set_tau_e(byref(self), c_double(value))
        @property
        def tau_inc(self):
            clibreboundx.rebx_get_tau_inc.restype = c_double
            return clibreboundx.rebx_get_tau_inc(byref(self))
        @tau_inc.setter
        def tau_inc(self, value):
            clibreboundx.rebx_set_tau_inc(byref(self), c_double(value))
        @property
        def tau_omega(self):
            clibreboundx.rebx_get_tau_omega.restype = c_double
            return clibreboundx.rebx_get_tau_omega(byref(self))
        @tau_omega.setter
        def tau_omega(self, value):
            clibreboundx.rebx_set_tau_omega(byref(self), c_double(value))
        @property
        def tau_Omega(self):
            clibreboundx.rebx_get_tau_Omega.restype = c_double
            return clibreboundx.rebx_get_tau_Omega(byref(self))
        @tau_Omega.setter
        def tau_Omega(self, value):
            clibreboundx.rebx_set_tau_Omega(byref(self), c_double(value))

        rebound.Particle.tau_a = tau_a
        rebound.Particle.tau_e = tau_e
        rebound.Particle.tau_inc = tau_inc
        rebound.Particle.tau_omega = tau_omega
        rebound.Particle.tau_Omega = tau_Omega

# Need to put fields after class definition because of self-referencing
Extras._fields_ = [("sim", POINTER(Simulation)),
                ("forces", POINTER(CFUNCTYPE(None, POINTER(Simulation)))),
                ("ptm", POINTER(CFUNCTYPE(None, POINTER(Simulation)))),
                ("Nforces", c_int),
                ("Nptm", c_int),
                ("modify_orbits_forces", rebx_params_modify_orbits_forces),
                ("modify_orbits_direct", rebx_params_modify_orbits_direct),
                ("gr", rebx_params_gr)]

