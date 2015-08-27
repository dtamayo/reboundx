from ctypes import *
import os
import rebound

#Find the migration C library
pymodulespath = os.path.dirname(__file__)
try:
    clibreboundxf = CDLL(pymodulespath + '/../libreboundxf.so', RTLD_GLOBAL)
except:
    print("Cannot find library 'libreboundxf.so'.")
    raise

class rebxf_param_elements_forces(Structure):
    _fields_ = [("allocatedN", c_int),
                ("tau_a", POINTER(c_double)),
                ("tau_e", POINTER(c_double)),
                ("tau_inc", POINTER(c_double)),
                ("tau_omega", POINTER(c_double)),
                ("e_damping_p", c_double)]                     

class rebxf_param_elements_direct(Structure):
    _fields_ = [("allocatedN", c_int),
                ("tau_a", POINTER(c_double)),
                ("tau_e", POINTER(c_double)),
                ("tau_inc", POINTER(c_double)),
                ("tau_omega", POINTER(c_double)),
                ("e_damping_p", c_double)]                     

class rebxf_param_gr(Structure):
    _fields_ = [("c", c_double)]

class rebxf_params(Structure):
    _fields_ = [("sim", reb_simulation)

def mod_test():
    return clibreboundxf.rebxf_modify_elements
def modify_elements():
    func_address = cast(clibreboundxf.rebxf_modify_elements, c_void_p).value
    return func_address

def forces():
    func_address = cast(clibreboundxf.rebxf_forces, c_void_p).value
    return func_address

def addxf(sim):
    clibreboundxf.rebxf_addxf(sim.simulation)

class Params(object):
    simulation = None
    params = None

    def __init__(self, sim):
        self.simulation = sim.simulation
        clibreboundxf.rebxf_addxf.restype = POINTER(rebxf_params)
        self.params = clibreboundxf.rebxf_addxf(self.simulation)

    #TODO: find a way to set individual elements from python, e.g., xf.tau_a[2] = 1.e3
    def __del__(self):
        clibreboundxf.rebxf_free_xfparams(self.params)

    @property
    def tau_a(self):
        return [self.params.contents.tau_a[i] for i in range(self.simulation.contents.N)]
    @tau_a.setter
    def tau_a(self, tau_a):
        if len(tau_a) is not self.simulation.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_a))(*tau_a)
        clibreboundxf.rebxf_set_tau_a(self.simulation, byref(arr))

    @property
    def tau_e(self):
        return [self.params.contents.tau_e[i] for i in range(self.simulation.contents.N)]
    @tau_e.setter
    def tau_e(self, tau_e):
        if len(tau_e) is not self.simulation.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_e))(*tau_e)
        clibreboundxf.rebxf_set_tau_e(self.simulation, byref(arr))

    @property
    def tau_inc(self):
        return [self.params.contents.tau_inc[i] for i in range(self.simulation.contents.N)]
    @tau_inc.setter
    def tau_inc(self, tau_inc):
        if len(tau_inc) is not self.simulation.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_inc))(*tau_inc)
        clibreboundxf.rebxf_set_tau_inc(self.simulation, byref(arr))

    @property
    def tau_pomega(self):
        return [self.params.contents.tau_pomega[i] for i in range(self.simulation.contents.N)]
    @tau_pomega.setter
    def tau_pomega(self, tau_pomega):
        if len(tau_pomega) is not self.simulation.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_pomega))(*tau_pomega)
        clibreboundxf.rebxf_set_tau_pomega(self.simulation, byref(arr))

'''def set_e_damping_p(value):
    clibreboundxf.rebxf_set_e_damping_p(c_double(value))

def set_migration(tau_a):
    arr = (c_double * len(tau_a))(*tau_a)
    clibreboundxf.rebxf_set_migration(byref(arr),c_int(len(tau_a)))

def set_e_damping(tau_e, p=0.):
    c_double.in_dll(libreboundxf, "e_damping_p").value = p
    arr = (c_double * len(tau_e))(*tau_e)
    clibreboundxf.rebxf_set_e_damping(byref(arr),c_int(len(tau_e)))
    
def set_e_damping(tau_e):
    arr = (c_double * len(tau_e))(*tau_e)
    clibreboundxf.rebxf_set_e_damping(byref(arr),c_int(len(tau_e)))
    
def set_i_damping(tau_i):
    arr = (c_double * len(tau_i))(*tau_i)
    clibreboundxf.rebxf_set_i_damping(byref(arr),c_int(len(tau_i)))

def set_peri_precession(tau_po):
    arr = (c_double * len(tau_po))(*tau_po)
    clibreboundxf.rebxf_set_peri_precession(byref(arr), c_int(len(tau_po)))
   
def set_peri_precession(gamma, Rc, podot):
    clibreboundxf.set_peri_precession(c_double(gamma), c_double(Rc), c_double(podot))
'''

