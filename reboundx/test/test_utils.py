from ctypes import c_void_p, c_double
from reboundx import clibreboundx


clibreboundx.rebx_geopotential_potential.argtypes = [c_void_p, c_double, c_double, c_double]
clibreboundx.rebx_geopotential_potential.restype = c_double


def geopotential_value(model, phi, r, lambda_body):
    """ Calls directly to rebx_geopotential_potential via ctypes, 
    whithout passing by the forces loop of REBOUND. It is useful to 
    evaluate the model in a random point, isolated from the rest of 
    the simulation. """
    return clibreboundx.rebx_geopotential_potential(
        model._ptr, c_double(phi), c_double(r), c_double(lambda_body)
    )