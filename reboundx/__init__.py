# -*- coding: utf-8 -*-
#Make changes for python 2 and 3 compatibility
try:
    import builtins     # if this success it's python 3.x
    builtins.xrange = range
    builtins.basestring = (str,bytes)
except ImportError:
    pass                # python 2.x

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

# Get string for build in libreboundx.so
def build_date():
    return c_char_p.in_dll(clibreboundx, "rebx_build_str").value.decode("ascii")
def build_version():
    return c_char_p.in_dll(clibreboundx, "rebx_version_str").value.decode("ascii")
#Check for version
try:
    moduleversion = pkg_resources.require("reboundx")[0].version
    libreboundxversion = build_version()
    if moduleversion != libreboundxversion:
        print("WARNING: python module and libreboundx have different version numbers: '%s' vs '%s'.\n".format(moduleversion, libreboundxversion))
except:
    pass    # this check fails in python 3. Problem with setuptools

# When reboundx is imported, first monkey patch rebound.Particle so that we can add new parameters to the particles:

import rebound
def monkeyset(self, name, value):
    if (name not in rebound.Particle.__dict__) and (not hasattr(super(rebound.Particle, self), name)):
        # create new param in c
        clibreboundx.rebx_set_param_double(byref(self), c_char_p(name.encode('utf-8')), c_double(value))
    else:
        rebound.Particle.default_set(self, name, value)
        
def monkeyget(self, name):
    if (name not in rebound.Particle.__dict__) and (not hasattr(super(rebound.Particle, self), name)):
        # check param in c
        clibreboundx.rebx_get_param_double.restype = c_double
        return clibreboundx.rebx_get_param_double(byref(self), c_char_p(name.encode('utf-8')))
    else:
        return super(rebound.Particle, self).__getattr__(name)

rebound.Particle.default_set = rebound.Particle.__setattr__
rebound.Particle.__setattr__ = monkeyset
rebound.Particle.__getattr__ = monkeyget

from .extras import Extras
from .tools import coordinates, install_test
__all__ = ["Extras", "coordinates", "install_test"]
