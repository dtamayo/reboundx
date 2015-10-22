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
def build_str():
    return c_charp_p.in_dll(clibreboundx, "rebx_build_str").value.decode("ascii")

#Check for version
try:
    moduleversion = pkg_resources.require("reboundx")[0].version
    libreboundxversion = build_str()
    if moduleversion != libreboundxversion:
        print("WARNING: python module and libreboundx have different version numbers: '%s' vs '%s'.\n".format(moduleversion, libreboundxversion))
except:
    pass    # this check fails in python 3. Problem with setuptools

from .extras import Extras
__all__ = ["Extras"]
