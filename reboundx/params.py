import rebound
import sys
if sys.version_info[:2] >= (3, 8):
    from collections.abc import MutableMapping
else:
    from collections import MutableMapping
from .extras import Param, Node, Force, Operator, Extras, REBX_CTYPES
from . import clibreboundx
from ctypes import byref, c_double, c_int, c_int32, c_int64, c_uint, c_uint32, c_longlong, c_char_p, POINTER, cast
from ctypes import c_void_p, memmove, sizeof, addressof
from rebound.tools import hash as rebhash

class Params(MutableMapping):
    def __init__(self, parent):
        self.verbose = 0        # set to 1 to diagnose problems
        self.parent = parent    # Particle, Force, Operator. Will work with any ctypes.Structure with appropriate ._sim and .ap fields

        offset = type(parent).ap.offset # Need this hack to get address of initially NULL ap ptr. See my stackoverflow
        self.ap = (c_void_p).from_buffer(parent, offset)
      
        # We want to be able to access params from objects like particles, objects etc. These must somehow point back to rebx for memory management
        # We can't add a rebx pointer in rebound.Particle, but can use its sim pointer. So we add a sim pointer to all rebx objects that need params.
        extrasvp = parent._sim.contents.extras
        if not extrasvp: # .extras = None
            raise AttributeError("Need to attach reboundx.Extras instance to simulation before setting params.")
        else:
            self.rebx = cast(extrasvp, POINTER(Extras))

    def __getitem__(self, key):
        param_type = clibreboundx.rebx_get_type(self.rebx, c_char_p(key.encode('ascii')))
        ctype = REBX_CTYPES[param_type]
        if ctype == None:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' not found in REBOUNDx. Need to register it first.".format(key))
        
        clibreboundx.rebx_get_param.restype = c_void_p
        valptr = clibreboundx.rebx_get_param(self.rebx, self.ap, c_char_p(key.encode('ascii')))

        if ctype == c_void_p: # Don't know how to cast it, so return for user to cast
            if valptr is None:
                raise AttributeError("REBOUNDx Error: Parameter '{0}' not found on object.".format(key))
            return valptr
        
        valptr = cast(valptr, POINTER(ctype))
        try:
            val = valptr.contents.value # return python int or float rather than c_int or c_double
        except AttributeError:
            val = valptr.contents # Structure, return ctypes object
        except ValueError: # NULL access
            raise AttributeError("REBOUNDx Error: Parameter '{0}' not found on object.".format(key))

        return val

    def __setitem__(self, key, value):
        param_type = clibreboundx.rebx_get_type(self.rebx, c_char_p(key.encode('ascii')))
        ctype = REBX_CTYPES[param_type]
        if ctype == None:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' not found in REBOUNDx. Need to register it first.".format(key))
        if ctype == c_double:
            clibreboundx.rebx_set_param_double(self.rebx, byref(self.ap), c_char_p(key.encode('ascii')), c_double(value))
        if ctype == c_int:
            clibreboundx.rebx_set_param_int(self.rebx, byref(self.ap), c_char_p(key.encode('ascii')), c_int(value))
        if ctype == c_uint32:
            clibreboundx.rebx_set_param_uint32(self.rebx, byref(self.ap), c_char_p(key.encode('ascii')), value)

        if ctype == Force:
            if not isinstance(value, Force):
                raise AttributeError("REBOUNDx Error: Parameter '{0}' must be assigned a Force object.".format(key))
            clibreboundx.rebx_set_param_pointer(self.rebx, byref(self.ap), c_char_p(key.encode('ascii')), byref(value))

        if ctype == rebound.Orbit:
            if not isinstance(value, rebound.Orbit):
                raise AttributeError("REBOUNDx Error: Parameter '{0}' must be assigned an Orbit object.".format(key))
            clibreboundx.rebx_set_param_pointer(self.rebx, byref(self.ap), c_char_p(key.encode('ascii')), byref(value))
        if ctype == c_void_p:
            clibreboundx.rebx_set_param_pointer(self.rebx, byref(self.ap), c_char_p(key.encode('ascii')), byref(value))

    def __delitem__(self, key):
        raise AttributeError("REBOUNDx Error: Removing particle params not implemented.")

    def __iter__(self):
        raise AttributeError("REBOUNDx Error: Iterator for params not implemented.")

    def __len__(self):
        clibreboundx.rebx_len.restype = c_int
        return clibreboundx.rebx_len(self.ap)
