# -*- coding: utf-8 -*-
import warnings

#Find suffix
import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

#Import shared C library
import os
import sys
pymodulespath = os.path.dirname(os.path.abspath(__file__))
pymodulespath = os.path.abspath(os.path.join(pymodulespath, os.pardir))
__libpath__ = os.path.join(pymodulespath, "libreboundx"+suffix)

# On Windows, libreboundx.pyd depends on librebound.pyd. Unlike Unix's
# RPATH/RUNPATH, Windows resolves dependent DLLs via the caller's working
# directory, PATH, or explicitly-registered "DLL directories". rebound is
# installed in site-packages, which by default isn't on the DLL search
# path, so we register it via os.add_dll_directory() before the load.
if sys.platform == 'win32' and hasattr(os, 'add_dll_directory'):
    import sysconfig as _sc
    _sitepkg = _sc.get_path('platlib')
    if _sitepkg and os.path.isdir(_sitepkg):
        os.add_dll_directory(_sitepkg)

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
