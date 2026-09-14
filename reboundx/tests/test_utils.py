from ctypes import c_void_p, c_double, POINTER, byref
from reboundx import clibreboundx


clibreboundx.rebx_spherical_harmonics_potential.argtypes = [c_void_p, c_double, c_double, c_double]
clibreboundx.rebx_spherical_harmonics_potential.restype = c_double


def spherical_harmonics_value(model, phi, r, lambda_body):
    """ Calls directly to rebx_spherical_harmonics_potential via ctypes, 
    whithout passing by the forces loop of REBOUND. It is useful to 
    evaluate the model in a random point, isolated from the rest of 
    the simulation. """
    return clibreboundx.rebx_spherical_harmonics_potential(
        model._ptr, c_double(phi), c_double(r), c_double(lambda_body)
    )

clibreboundx.rebx_spherical_harmonics_acceleration.argtypes = [
    c_void_p, c_double, c_double, c_double,
    POINTER(c_double), POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_double),
]
clibreboundx.rebx_spherical_harmonics_acceleration.restype = None


def spherical_harmonics_acceleration(model, phi, r, lambda_body):
    """ Calls directly to rebx_spherical_harmonics_acceleration via ctypes,
    isolated from the REBOUND integrator/timestep. Returns (a_r, a_phi, a_lambda)
    in the body-fixed spherical frame. """
    a_r = c_double()
    a_phi = c_double()
    a_lambda = c_double()
    sinphi = c_double()
    cosphi = c_double()

    clibreboundx.rebx_spherical_harmonics_acceleration(
        model._ptr, c_double(phi), c_double(r), c_double(lambda_body),
        byref(a_r), byref(a_phi), byref(a_lambda),
        byref(sinphi), byref(cosphi),
    )
    return a_r.value, a_phi.value, a_lambda.value