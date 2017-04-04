# -*- coding: utf-8 -*-

#Find suffix
import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

#Import shared C library
import os
pymodulespath = os.path.dirname(__file__)
from ctypes import *
clibreboundx = cdll.LoadLibrary(pymodulespath + '/../libreboundx' + suffix)

# Version
__version__ = c_char_p.in_dll(clibreboundx, "rebx_version_str").value.decode('ascii')

# Build
__build__ = c_char_p.in_dll(clibreboundx, "rebx_build_str").value.decode('ascii')
# Check for version

# Githash
__githash__ = c_char_p.in_dll(clibreboundx, "rebx_githash_str").value.decode('ascii')

try:
    moduleversion = pkg_resources.require("reboundx")[0].version
    libreboundxversion = __version__
    if moduleversion != libreboundxversion:
        print("WARNING: python module and libreboundx have different version numbers: '%s' vs '%s'.\n".format(moduleversion, libreboundxversion))
except:
    pass    # this check fails in python 3. Problem with setuptools

# Monkeypatch rebound.Particle to have a params property
import rebound

@property
def params(self):
    params = Params(self)
    return params

rebound.Particle.params = params

from .extras import Extras, Param, Effect
from .tools import coordinates, install_test
from .params import Params
__all__ = ["__version__", "__build__", "__githash__", "Extras", "Param", "Effect", "Params", "coordinates", "install_test"]
