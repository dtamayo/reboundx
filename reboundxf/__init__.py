from ctypes import *
import os
import sys
sys.path.append('../../rebound/rebound/')
import rebound

#Find the migration C library
pymodulespath = os.path.dirname(__file__)
try:
    clibreboundxf = CDLL(pymodulespath + '/../libreboundxf.so', RTLD_GLOBAL)
except:
    print("Cannot find library 'libreboundxf.so'.")
    raise

def modify_elements():
    return clibreboundxf.rebxf_modify_elements

def forces():
    return clibreboundxf.rebxf_forces

def addxf(sim):
    clibreboundxf.rebxf_addxf(sim.simulation)

class Params(object):
    simulation = None
    def __init__(self, sim):
        self.simulation = sim.simulation
        clibreboundxf.rebxf_addxf(self.simulation)

    @property
    def tau_a(self):
        clibreboundxf.rebxf_get_tau_a.restype = POINTER(c_double)
        tau_a = clibreboundxf.rebxf_get_tau_a(self.simulation) 
        return [tau_a[i] for i in range(self.simulation.contents.N)]
    @tau_a.setter
    def tau_a(self, tau_a):
        if len(tau_a) is not self.simulation.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_a))(*tau_a)
        clibreboundxf.rebxf_set_tau_a(self.simulation, byref(arr))

    @property
    def tau_e(self):
        clibreboundxf.rebxf_get_tau_e.restype = POINTER(c_double)
        tau_e = clibreboundxf.rebxf_get_tau_e(self.simulation) 
        return [tau_e[i] for i in range(self.simulation.contents.N)]
    @tau_e.setter
    def tau_e(self, tau_e):
        if len(tau_e) is not self.simulation.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_e))(*tau_e)
        clibreboundxf.rebxf_set_tau_e(self.simulation, byref(arr))

    @property
    def tau_inc(self):
        clibreboundxf.rebxf_get_tau_inc.restype = POINTER(c_double)
        tau_inc = clibreboundxf.rebxf_get_tau_inc(self.simulation) 
        return [tau_inc[i] for i in range(self.simulation.contents.N)]
    @tau_inc.setter
    def tau_inc(self, tau_inc):
        if len(tau_inc) is not self.simulation.contents.N:
            raise AttributeError('You must pass an array of timescales with as many elements as there are particles in the simulation')
        arr = (c_double * len(tau_inc))(*tau_inc)
        clibreboundxf.rebxf_set_tau_inc(self.simulation, byref(arr))

    @property
    def tau_pomega(self):
        clibreboundxf.rebxf_get_tau_pomega.restype = POINTER(c_double)
        tau_pomega = clibreboundxf.rebxf_get_tau_pomega(self.simulation) 
        return [tau_pomega[i] for i in range(self.simulation.contents.N)]
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

