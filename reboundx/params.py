import rebound
from collections import MutableMapping
from .extras import rebx_param
from . import clibreboundx
from ctypes import byref, c_double, c_int, c_char_p, POINTER, cast

class Params(MutableMapping):
    def __init__(self, parent):
        # First make sure a rebx instance is attached to sim, otherwise memory isn't allocated.
        if parent.__class__ == rebound.particle.Particle: 
            if not parent._sim.contents.extras:
                raise AttributeError("Need to attach reboundx.Extras instance to simulation before setting particle params.")

        self.parent = parent
        self.type_names = {"float":"double", "int":"int"}
        self.c_types = {"double":c_double, "int":c_int}

    def __getitem__(self, key):
        name = rebound.hash(key).value
        current = cast(self.parent.ap, POINTER(rebx_param))
        while current:
            param = current.contents
            if param.hash == name:
                res = None
                for data_type in self.c_types.keys():
                    if param.type_hash == rebound.hash(data_type).value:
                        res = cast(param.paramPtr, POINTER(self.c_types[data_type]))
                if res is None:
                    raise AttributeError("Data type for {0} param not supported.  Need to update reboundx/params.py".format(key)) 
                else:
                    return res.contents.value
            current = param.next

        raise AttributeError("{0} parameter not found.".format(key))
                
    def __setitem__(self, key, value):
        type_name = self.type_names[type(value).__name__]
        c_type = self.c_types[type_name]
        function = getattr(clibreboundx, "rebx_set_param_"+type_name)
        function(byref(self.parent), c_char_p(key.encode('ascii')), c_type(value))
   
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
