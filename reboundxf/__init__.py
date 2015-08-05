from ctypes import *
import os

#Find the migration C library
pymodulespath = os.path.dirname(__file__)
try:
    clibreboundxf = CDLL(pymodulespath + '/../libreboundxf.so', RTLD_GLOBAL)
except:
    print("Cannot find library 'libreboundxf.so'.")
    raise

def modify_elements():
    return clibreboundxf.reboundxf_modify_elements

def forces():
    return clibreboundxf.reboundxf_forces

def set_e_damping_p(value):
    clibreboundxf.reboundxf_set_e_damping_p(c_double(value))

def set_migration(tau_a):
    arr = (c_double * len(tau_a))(*tau_a)
    clibreboundxf.reboundxf_set_migration(byref(arr),c_int(len(tau_a)))

def set_e_damping(tau_e, p=0.):
    c_double.in_dll(libreboundxf, "e_damping_p").value = p
    arr = (c_double * len(tau_e))(*tau_e)
    clibreboundxf.reboundxf_set_e_damping(byref(arr),c_int(len(tau_e)))
    
def set_e_damping(tau_e):
    arr = (c_double * len(tau_e))(*tau_e)
    clibreboundxf.reboundxf_set_e_damping(byref(arr),c_int(len(tau_e)))
    
def set_i_damping(tau_i):
    arr = (c_double * len(tau_i))(*tau_i)
    clibreboundxf.reboundxf_set_i_damping(byref(arr),c_int(len(tau_i)))

def set_peri_precession(tau_po):
    arr = (c_double * len(tau_po))(*tau_po)
    clibreboundxf.reboundxf_set_peri_precession(byref(arr), c_int(len(tau_po)))

'''    
def set_peri_precession(gamma, Rc, podot):
    clibreboundxf.set_peri_precession(c_double(gamma), c_double(Rc), c_double(podot))
'''

