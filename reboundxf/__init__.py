from ctypes import *
import os
'''
#Find the migration C library
pymodulespath = os.path.dirname(__file__)
try:
    clibreboundxf = CDLL(pymodulespath + '/../libreboundxf.so', RTLD_GLOBAL)
except:
    print("Cannot find library 'libreboundxf.so'.")
    raise
'''
try:
    import builtins      # if this succeeds it's python 3.x
    builtins.xrange = range
    builtins.basestring = (str,bytes)
except ImportError:
    pass                 # python 2.x

def forces():
    return clibreboundxf.disk_forces

def init(N):
    clibreboundxf.init(c_int(N))
    

def set_migration(tau_a):
    arr = (c_double * len(tau_a))(*tau_a)
    clibreboundxf.set_migration(byref(arr))

def set_e_damping(tau_e, p=0.):
    c_double.in_dll(libreboundxf, "e_damping_p").value = p
    arr = (c_double * len(tau_e))(*tau_e)
    clibreboundxf.set_e_damping(byref(arr))
    
def set_e_damping(tau_e):
    arr = (c_double * len(tau_e))(*tau_e)
    clibreboundxf.set_e_damping(byref(arr))
    
def set_i_damping(tau_i):
    arr = (c_double * len(tau_i))(*tau_i)
    clibreboundxf.set_i_damping(byref(arr))
    
def set_peri_precession(gamma, Rc, podot):
    clibreboundxf.set_peri_precession(c_double(gamma), c_double(Rc), c_double(podot))

