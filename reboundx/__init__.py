# -*- coding: utf-8 -*-
import warnings

#Find suffix
import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

#Import shared C library
import os
pymodulespath = os.path.dirname(os.path.abspath(__file__))
pymodulespath = os.path.abspath(os.path.join(pymodulespath, os.pardir))
__libpath__ = os.path.join(pymodulespath, "libreboundx"+suffix)
from ctypes import *
clibreboundx = cdll.LoadLibrary(__libpath__)

# Version
__version__ = c_char_p.in_dll(clibreboundx, "rebx_version_str").value.decode('ascii')

# Build
__build__ = c_char_p.in_dll(clibreboundx, "rebx_build_str").value.decode('ascii')
# Check for version

# Githash
__githash__ = c_char_p.in_dll(clibreboundx, "rebx_githash_str").value.decode('ascii')

try:
    import pkg_resources
    moduleversion = pkg_resources.require("reboundx")[0].version
    libreboundxversion = __version__
    if moduleversion != libreboundxversion:
        warnings.warn("WARNING: python module and libreboundx have different version numbers: '%s' vs '%s'.\n" % (moduleversion, libreboundxversion), ImportWarning)
except:
    pass    # this check fails in python 3. Problem with setuptools

# Monkeypatch rebound.Particle to have a params property
import rebound

@property
def params(self):
    params = Params(self)
    return params

rebound.Particle.params = params

from .extras import Extras, Param, Node, Force, Operator, integrators, Interpolator
from .simulationarchive import Simulationarchive
from .tools import coordinates, install_test
from .params import Params

__all__ = ["__version__", "__build__", "__githash__", "Extras", "Simulationarchive", "Param", "Interpolator", "Params", "coordinates", "integrators"]

from ctypes import c_int, c_double, c_void_p, c_char_p, POINTER
 
clibreboundx.rebx_spherical_harmonics_create.restype = c_void_p
clibreboundx.rebx_spherical_harmonics_create.argtypes = [
    POINTER(c_double),
    POINTER(c_double),
    c_double,
    c_double,
]
 
clibreboundx.rebx_spherical_harmonics_free.restype = None
clibreboundx.rebx_spherical_harmonics_free.argtypes = [c_void_p]
 
# Extras ya debería estar definida en este mismo __init__.py o en extras.py;
# si está en extras.py, mueve estas 3 líneas después de "from .extras import Extras"
clibreboundx.rebx_set_param_pointer.restype = None
clibreboundx.rebx_set_param_pointer.argtypes = [
    POINTER(Extras),
    c_void_p,   # dirección de particle.ap (obtenida con byref(particle, Particle.ap.offset))
    c_char_p,
    c_void_p,
]
 
 

from .spherical_harmonics import spherical_harmonics_model
__all__ += ["spherical_harmonics_model"]