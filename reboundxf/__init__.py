from ctypes import *
import os

#Find the migration C library
pymodulespath = os.path.dirname(__file__)
try:
    clibreboundxf = CDLL(pymodulespath + '/../reboundxf.so', RTLD_GLOBAL)
except:
    print("Cannot find library 'libreboundxf.so'.")
    raise

def forces():
    return clibreboundxf.forces

def set_migration(tau_a):
    arr = (c_double * len(tau_a))(*tau_a)
    clibreboundxf.set_migration(byref(arr),c_int(len(tau_a)))

def set_e_damping(tau_e, p=0.):
    c_double.in_dll(libreboundxf, "e_damping_p").value = p
    arr = (c_double * len(tau_e))(*tau_e)
    clibreboundxf.set_e_damping(byref(arr),c_int(len(tau_e)))
    
def set_e_damping(tau_e):
    arr = (c_double * len(tau_e))(*tau_e)
    clibreboundxf.set_e_damping(byref(arr),c_int(len(tau_e)))
    
def set_i_damping(tau_i):
    arr = (c_double * len(tau_i))(*tau_i)
    clibreboundxf.set_i_damping(byref(arr),c_int(len(tau_i)))
    
'''not yet implemented
def set_peri_precession(gamma, Rc, podot):
    clibreboundxf.set_peri_precession(c_double(gamma), c_double(Rc), c_double(podot))
'''

