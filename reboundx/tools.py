from . import clibreboundx
from ctypes import c_double

coordinates = {"JACOBI":0, "BARYCENTRIC":1, "HELIOCENTRIC":2} # to use C version's REBX_COORDINATES enum

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


