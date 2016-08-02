import rebound
from collections import MutableMapping
from .extras import rebx_param
from . import clibreboundx
from ctypes import byref, c_double, c_int, c_char_p, POINTER, cast
from numpy.ctypeslib import as_ctypes

class Params(MutableMapping):
    def __init__(self, parent):
        # First make sure a rebx instance is attached to sim, otherwise memory isn't allocated.
        if parent.__class__ == rebound.particle.Particle: 
            if not parent._sim.contents.extras:
                raise AttributeError("Need to attach reboundx.Extras instance to simulation before setting particle params.")

        self.parent = parent

    def __getitem__(self, key):
        name = rebound.hash(key).value
        current = cast(self.parent.ap, POINTER(rebx_param))
        while current:
            param = current.contents
            if param.hash == name:
                res = None
                if param.type_hash == rebound.hash('int').value:
                    res = cast(param.paramPtr, POINTER(c_int))
                elif param.type_hash == rebound.hash('double').value:
                    res = cast(param.paramPtr, POINTER(c_double))
                
                if res is None:
                    raise AttributeError("Data type for {0} param not supported.".format(key)) 
                else:
                    return res.contents.value
            current = param.next

        raise AttributeError("{0} parameter not found.".format(key))
                
    def __setitem__(self, key, value):
        typ = type(value)
        if(typ.__module__ == 'numpy'):
            ctypevar = as_ctypes(value)
            ctype = type(ctypevar)
            if ctype == c_int:
                clibreboundx.rebx_set_param_int(byref(self.parent), c_char_p(key.encode('ascii')), ctypevar)
                return
            elif ctype == c_double:
                clibreboundx.rebx_set_param_double(byref(self.parent), c_char_p(key.encode('ascii')), ctypevar)
                return
        else:
            if typ == int:
                clibreboundx.rebx_set_param_int(byref(self.parent), c_char_p(key.encode('ascii')), c_int(value))
                return
            elif typ == float:
                clibreboundx.rebx_set_param_double(byref(self.parent), c_char_p(key.encode('ascii')), c_double(value))
                return
        raise AttributeError("Data type for {0} param not supported.".format(key))

    def __delitem__(self, key):
        success = clibreboundx.rebx_remove_param(byref(self.parent), c_char_p(key.encode('ascii')))
        if not success:
            raise AttributeError("{0} parameter to delete not found.".format(key))

    def __iter__(self):
        raise AttributeError("Iterator for params not implemented.")

    def __len__(self):
        current = cast(self.parent.ap, POINTER(rebx_param))
        ctr = 0
        while current:
            ctr += 1
            current = current.contents.next
        return ctr
