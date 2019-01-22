import rebound
from collections import MutableMapping
from .extras import Param, Node, Force, Operator, Extras
from . import clibreboundx
from ctypes import byref, c_double, c_int, c_int32, c_int64, c_uint, c_uint32, c_longlong, c_char_p, POINTER, cast
from ctypes import c_void_p, memmove, sizeof, addressof
import numpy as np
from rebound.tools import hash as rebhash

# TODO: Trim this down to what's needed now that we don't support arrays

#### Update here to add new types. Update only rebx_param_types with new C rebx_types
#### and self.types for new python types.
#### ADD ONLY TO THE END OF EACH LIST 
        

rebx_ctypes = { 0: None,
                1: c_double,
                2: c_int,
                3: c_void_p}

class Params(MutableMapping):
    def __init__(self, parent):
        self.verbose = 0        # set to 1 to diagnose problems
        self.parent = parent    # Particle, Force, Effect. Will work with any ctypes.Structure with appropriate ._sim and .ap fields

        offset = type(parent).ap.offset # Need this hack to get address of initially NULL ap ptr. See my stackoverflow
        self.ap = (c_void_p).from_buffer(parent, offset)
       
        extrasvp = parent._sim.contents.extras
        if not extrasvp: # .extras = None
            raise AttributeError("Need to attach reboundx.Extras instance to simulation before setting params.")
        else:
            self.rebx = cast(extrasvp, POINTER(Extras))

    def __getitem__(self, key):
        param_type = clibreboundx.rebx_get_type(c_char_p(key.encode('ascii')))
        if param_type == 0:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' not found in REBOUNDx. Need to register it first.".format(key))
        if param_type == 3:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' is a pointer. Need to manually set pointers with rebx.add_pointer_param.".format(key))
        clibreboundx.rebx_get_param.restype = POINTER(Param)
        paramptr = clibreboundx.rebx_get_param(self.rebx, self.ap, c_char_p(key.encode('ascii')))
        ctype = rebx_ctypes[param_type]
        valptr = cast(paramptr.contents.value, POINTER(ctype))
        return valptr.contents.value
        '''    
            param = paramptr.contents
        else:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' not found.".format(key))

        try:
            pythontype = rebx_ptype[param.python_type]
        except KeyError:
            pythontype = rebx_defaultptype[param.param_type] # Parameter was set in C, use default

        dtype, ctype, rebxtype = rebx_types[pythontype]
        if self.verbose:
            print("param.python_type= {0}, pythontype = {1}".format(param.python_type, pythontype))
            print("dtype = {0}, ctype = {1}, rebxtype = {2}".format(dtype, ctype, rebxtype))

        valptr = cast(param.value, POINTER(ctype))
        try: # TODO More explicit check for whether it's a Structure
            val = valptr.contents.value # cases where scalar is a ctype, so we get python type back 
            val = pythontype(val)
        except AttributeError:
            if self.verbose: print("Parameter is a Structure")
            val = valptr.contents # cases where scalar is a Structure, e.g., rebound.Orbit
    
        return val
        '''

    def __setitem__(self, key, value):
        param_type = clibreboundx.rebx_get_type(c_char_p(key.encode('ascii')))
        if param_type == 0:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' not found in REBOUNDx. Need to register it first.".format(key))
        if param_type == 3:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' is a pointer. Need to manually set pointers with rebx.add_pointer_param.".format(key))
        if param_type == 1:
            valptr = byref(c_double(value))
        if param_type == 2:
            valptr = byref(c_int(value))
        clibreboundx.rebx_set_param_copy(self.rebx, byref(self.ap), c_char_p(key.encode('ascii')), valptr)
        '''
        pythontype = type(value)
        try:
            dtype, ctype, rebxtype = rebx_types[pythontype]
            if self.verbose: print("ctype = {0}, rebxtype = {1}".format(ctype, rebxtype))
        except KeyError:
            raise AttributeError("REBOUNDx Error: Data type {0} for param '{1}' "
            "not supported.".format(pythontype, key))
       
        clibreboundx.rebx_get_param_struct.restype = POINTER(Param)
        paramptr = clibreboundx.rebx_get_param_struct(self.get_ap(), c_char_p(key.encode('ascii')))

        if paramptr:
            if self.verbose: print("GParameter {0} found.".format(key))
            if paramptr.contents.param_type != rebxtype:
                raise AttributeError("REBOUNDx Error: Can't update param '{0}' with "
                "new type {1}".format(key, pythontype))
        else: # parameter not found. Make new one
            if self.verbose: print("Parameter {0} not found. Making new one".format(key))
            clibreboundx.rebx_add_param(self.parent._sim, byref(self.get_ap()), c_char_p(key.encode('ascii')), rebxtype) 
            paramptr = clibreboundx.rebx_get_param_struct(self.get_ap(), c_char_p(key.encode('ascii')))
            if not paramptr:
                raise AttributeError("REBOUNDx Error: Couldn't add param '{0}'".format(key))
            paramptr.contents.python_type = c_int(rebx_penum[pythontype])

        val = cast(paramptr.contents.value, POINTER(ctype))
        val[0] = value
        '''
    def __delitem__(self, key):
        success = clibreboundx.rebx_remove_param(byref(self.get_ap()), c_char_p(key.encode('ascii')))
        if not success:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' to delete not found.".format(key))

    def __iter__(self):
        raise AttributeError("REBOUNDx Error: Iterator for params not implemented.")

    def __len__(self):
        clibreboundx.rebx_len.restype = c_int
        return clibreboundx.rebx_len(self.get_ap())

