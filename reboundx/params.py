import rebound
from collections import MutableMapping
from .extras import Param, Effect
from . import clibreboundx
from ctypes import byref, c_double, c_int, c_int32, c_int64, c_uint, c_uint32, c_longlong, c_char_p, POINTER, cast
from ctypes import c_void_p, memmove, sizeof, addressof
import numpy as np
from rebound.tools import hash as rebhash

# TODO: Trim this down to what's needed now that we don't support arrays

#### Update here to add new types. Update only rebx_param_types with new C rebx_types
#### and self.types for new python types.
#### ADD ONLY TO THE END OF EACH LIST 
        

rebx_param_types = [("REBX_TYPE_DOUBLE", float),
                    ("REBX_TYPE_INT", np.int32),
                    ("REBX_TYPE_UINT32", c_uint),
                    ("REBX_TYPE_ORBIT", rebound.Orbit),
                    ("REBX_TYPE_LONGLONG", np.int64),
                   ] # In same order as C enum, with default python types for each
val = {param_type[0]:value for (value, param_type) in enumerate(rebx_param_types)}

# tuples of corresponding (python type, dtype, ctype, rebxtype enum value). ONLY ADD TO THE END OF THIS LIST
# dtype should be same as python type for 'primitive' types, and object for any ctypes structure
type_tuples = [(int, int, c_int, val["REBX_TYPE_INT"]), 
               (np.int32, np.int32, c_int, val["REBX_TYPE_INT"]), 
               (float, float, c_double, val["REBX_TYPE_DOUBLE"]), 
               (np.float64, np.float64, c_double, val["REBX_TYPE_DOUBLE"]),
               (c_uint, object, c_uint, val["REBX_TYPE_UINT32"]),
               (c_uint32, object, c_uint32, val["REBX_TYPE_UINT32"]),
               (rebound.Orbit, object, rebound.Orbit, val["REBX_TYPE_ORBIT"]),
               (np.int64, np.int64, c_longlong, val["REBX_TYPE_LONGLONG"]),
              ]

python_types = [item[0] for item in type_tuples]
rebx_ptype = {enum:ptype for (enum, ptype) in enumerate(python_types)} # get python type from enum value
rebx_defaultptype = {val[rebxtype]:ptype for (rebxtype, ptype) in rebx_param_types} # get default python type from C rebx type
rebx_penum = {ptype:enum for (enum, ptype) in enumerate(python_types)} # get enum from python type
rebx_types = {item[0]:(item[1], item[2], item[3]) for item in type_tuples} # (dtype, ctype, rebxtype)from ptype

class Params(MutableMapping):
    def __init__(self, parent):
        self.verbose = 0 # set to 1 to diagnose problems
        self.parent = parent

        # Check rebx instance is attached, otherwise memory isn't allocated.
        if parent.__class__ == rebound.particle.Particle: 
            if not parent._sim.contents.extras:
                raise AttributeError("Need to attach reboundx.Extras instance to simulation before setting"
                        " particle params.")

    def __getitem__(self, key):
        clibreboundx.rebx_get_param_struct.restype = POINTER(Param)
        paramptr = clibreboundx.rebx_get_param_struct(self.get_ap(), c_char_p(key.encode('ascii')))
       
        if paramptr: 
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

    def __setitem__(self, key, value):
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
            if self.verbose: print("Parameter {0} found.".format(key))
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

    def __delitem__(self, key):
        success = clibreboundx.rebx_remove_param(byref(self.get_ap()), c_char_p(key.encode('ascii')))
        if not success:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' to delete not found.".format(key))

    def __iter__(self):
        raise AttributeError("REBOUNDx Error: Iterator for params not implemented.")

    def __len__(self):
        clibreboundx.rebx_len.restype = c_int
        return clibreboundx.rebx_len(self.get_ap())

    def get_ap(self):
        offset = type(self.parent).ap.offset
        return (c_void_p).from_buffer(self.parent, offset)


