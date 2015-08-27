from ctypes import *
import os
import rebound

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

class rebx_extras(Structure):
    _fields_ = [("sim", reb_simulation),
                ("forces", POINTER(CFUNCTYPE(None, POINTER(reb_simulation)))),
                ("ptm", POINTER(CFUNCTYPE(None, POINTER(reb_simulation)))),
                ("Nforces", c_int),
                ("Nptm", c_int),
                ("modify_orbits_forces", rebx_param_elements_forces),
                ("modify_orbits_direct", rebx_param_elements_direct),
                ("gr", rebx_param_gr)]

def mod_test():
    return clibreboundx.rebx_modify_elements
def modify_elements():
    func_address = cast(clibreboundx.rebx_modify_elements, c_void_p).value
    return func_address

def forces():
    func_address = cast(clibreboundx.rebx_forces, c_void_p).value
    return func_address

def addx(sim):
    clibreboundx.rebx_addx(sim.simulation)

class Params(object):
    simulation = None
    params = None

    def __init__(self, sim):
        self.simulation = sim.simulation
        clibreboundx.rebx_addx.restype = POINTER(rebx_params)
        self.params = clibreboundx.rebx_addx(self.simulation)

    #TODO: find a way to set individual elements from python, e.g., x.tau_a[2] = 1.e3
    def __del__(self):
        clibreboundx.rebx_free_xparams(self.params)

    @property
    def tau_a(self):
        return [self.params.contents.tau_a[i] for i in range(self.simulation.contents.N)]
    @tau_a.setter
    def tau_a(self, tau_a):
        if len(tau_a) is not self.simulation.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_a))(*tau_a)
        clibreboundx.rebx_set_tau_a(self.simulation, byref(arr))

    @property
    def tau_e(self):
        return [self.params.contents.tau_e[i] for i in range(self.simulation.contents.N)]
    @tau_e.setter
    def tau_e(self, tau_e):
        if len(tau_e) is not self.simulation.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_e))(*tau_e)
        clibreboundx.rebx_set_tau_e(self.simulation, byref(arr))

    @property
    def tau_inc(self):
        return [self.params.contents.tau_inc[i] for i in range(self.simulation.contents.N)]
    @tau_inc.setter
    def tau_inc(self, tau_inc):
        if len(tau_inc) is not self.simulation.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_inc))(*tau_inc)
        clibreboundx.rebx_set_tau_inc(self.simulation, byref(arr))

    @property
    def tau_pomega(self):
        return [self.params.contents.tau_pomega[i] for i in range(self.simulation.contents.N)]
    @tau_pomega.setter
    def tau_pomega(self, tau_pomega):
        if len(tau_pomega) is not self.simulation.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_pomega))(*tau_pomega)
        clibreboundx.rebx_set_tau_pomega(self.simulation, byref(arr))

'''def set_e_damping_p(value):
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

def set_peri_precession(tau_po):
    arr = (c_double * len(tau_po))(*tau_po)
    clibreboundx.rebx_set_peri_precession(byref(arr), c_int(len(tau_po)))
   
def set_peri_precession(gamma, Rc, podot):
    clibreboundx.set_peri_precession(c_double(gamma), c_double(Rc), c_double(podot))
'''

