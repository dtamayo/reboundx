import rebound
from collections import MutableMapping
from .extras import Param, Effect
from . import clibreboundx
from ctypes import byref, c_double, c_int, c_char_p, POINTER, cast, c_void_p, memmove, sizeof
from numpy.ctypeslib import as_ctypes
import numpy as np

class Params(MutableMapping):
    def __init__(self, parent):
        # First make sure a rebx instance is attached to sim, otherwise memory isn't allocated.
        if parent.__class__ == rebound.particle.Particle: 
            if not parent._sim.contents.extras:
                raise AttributeError("Need to attach reboundx.Extras instance to simulation before setting particle params.")

        self.parent = parent
        rebx_param_types = ["REBX_TYPE_DOUBLE", "REBX_TYPE_INT"]
        val = {param_type:value for (value, param_type) in enumerate(rebx_param_types)}

        self.param_types = {c_int:val["REBX_TYPE_INT"], c_double:val["REBX_TYPE_DOUBLE"]}   # ctype to rebx type enum value
        self.ctypes = {value:key for (key, value) in self.param_types.items()}              # rebx type enum value to ctype
        self.ctypesP = {int:c_int, float:c_double}                                          # python type to ctype
        self.python_types = {val["REBX_TYPE_INT"]:int, val["REBX_TYPE_DOUBLE"]:float}       # param_type to python type

    def as_ctypes(self, key, value):
        try:
            return self.ctypesP[type(value)](value)
        except KeyError:
            raise AttributeError("REBOUNDx Error: Data type {0} for param '{1}' not supported.".format(type(value), key))

    def __getitem__(self, key):
        name = rebound.hash(key).value
        clibreboundx.rebx_get_param_node.restype = POINTER(Param)
        node = clibreboundx.rebx_get_param_node(byref(self.parent), c_char_p(key.encode('ascii')))
        if node:
            param = node.contents
            if param.ndim == 1 and param.shape[0] == 1:
                data = cast(param.contents, POINTER(self.ctypes[param.param_type]))
                return data.contents.value
            else:
                ArrayType = self.ctypes[param.param_type]*param.size
                data = cast(param.contents, POINTER(ArrayType))
                array = np.frombuffer(data.contents, dtype=self.python_types[param.param_type], count=param.size)
                if param.ndim > 1:
                    shape = (param.shape[i] for i in range(param.ndim))
                    array = np.reshape(array, tuple(shape))
                return array
        else:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' not found.".format(key))

    def __setitem__(self, key, value):
        if type(value) == np.ndarray:
            clibreboundx.rebx_add_param_.restype = c_void_p
            val = clibreboundx.rebx_add_param_(byref(self.parent), c_char_p(key.encode('ascii')), 0, value.ndim, value.ctypes.shape_as(c_int))
            val = cast(val, POINTER(c_double))
            memmove(val, value.ctypes.data_as(POINTER(c_double)), sizeof(c_double)*value.size)
            return
        if(type(value).__module__ == 'numpy'):
            ctypevar = as_ctypes(value)
        else:
            ctypevar = self.as_ctypes(key, value)
        ctype = type(ctypevar)

        clibreboundx.rebx_get_param_node.restype = POINTER(Param)
        nodeptr = clibreboundx.rebx_get_param_node(byref(self.parent), c_char_p(key.encode('ascii')))

        if nodeptr:
            node = nodeptr.contents 
            if node.param_type != self.param_types[ctype]:
                raise AttributeError("REBOUNDx Error: Cannot update param '{0}' with incompatible data type '{1}'".format(key, ctype))
            val = node.contents
        else: # parameter not found. Make new one
            try:
                param_type = self.param_types[ctype]
            except KeyError:
                raise AttributeError("REBOUNDx Error: Data type {0} for param '{1}' not supported.".format(type(value), key))

            clibreboundx.rebx_add_param.restype = c_void_p
            val = clibreboundx.rebx_add_param(byref(self.parent), c_char_p(key.encode('ascii')), param_type, 1)

        val = cast(val, POINTER(ctype))
        val[0] = value
        return

    def __delitem__(self, key):
        success = clibreboundx.rebx_remove_param(byref(self.parent), c_char_p(key.encode('ascii')))
        if not success:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' to delete not found.".format(key))

    def __iter__(self):
        raise AttributeError("REBOUNDx Error: Iterator for params not implemented.")

    def __len__(self):
        current = cast(self.parent.ap, POINTER(Param))
        ctr = 0
        while current:
            ctr += 1
            current = current.contents.next
        return ctr
