from . import clibreboundx
from ctypes import c_double

coordinates = {"JACOBI":0, "BARYCENTRIC":1, "HELIOCENTRIC":2} # to use C version's REBX_COORDINATES enum

#make sure c is passed when it needs to be
def check_c(rebx, c):
    if c is not None: # user passed c explicitly
        return c
  
    # c was not passed by user
     
    if rebx.sim.contents.G == 1: # if G = 1 (default) return default c
        return 10064.915 # speed of light in AU, yr/2pi
    else:
        raise ValueError("If you change G, you must pass c (speed of light) in appropriate units to add_gr, add_gr_potential, add_gr_full, and radiation_forces.  Setting the units in the simulation does not work with REBOUNDx.  See ipython_examples/GeneralRelativity.ipynb and ipython_examples/Radiation_Forces_Debris_Disk.ipynb")

#function to test whether REBOUND shared library can be located and called correctly
def install_test():
    e = None
    try:
        clibreboundx.install_test.restype = c_double
    except Exception as e:
        return e
    try:
        x = clibreboundx.install_test()
    except Exception as e:
        return e
    try:
        if abs(x-0.17599665767) > 1.e-6:
            return 'Did not integrate to correct value of particles[1].x'
    except Exception as e:
        return e
    
    return e


