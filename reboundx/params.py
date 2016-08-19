import rebound
from collections import MutableMapping
from .extras import Param, Effect
from . import clibreboundx
from ctypes import byref, c_double, c_int, c_int32, c_int64, c_uint, c_uint32, c_char_p, POINTER, cast, c_void_p, memmove, sizeof
import numpy as np

class Params(MutableMapping):
    def __init__(self, parent):
        # First make sure a rebx instance is attached to sim, otherwise memory isn't allocated.
        if parent.__class__ == rebound.particle.Particle: 
            if not parent._sim.contents.extras:
                raise AttributeError("Need to attach reboundx.Extras instance to simulation before setting particle params.")

        self.parent = parent
        
        #### Update here to add new types #####
        
        rebx_param_types = ["REBX_TYPE_DOUBLE", "REBX_TYPE_INT", "REBX_TYPE_UINT32"]
        val = {param_type:value for (value, param_type) in enumerate(rebx_param_types)}

        # takes rebxtype and gives ctype and dtype for numpy arrays
        self.from_type =    {   val["REBX_TYPE_INT"]:(c_int, 'int32'), # REBOUNDx uses 32 bit int for all integer parameters
                                val["REBX_TYPE_DOUBLE"]:(c_double, 'float'),
                                val["REBX_TYPE_UINT32"]:(c_uint32, 'uint32'),
                            }
        # takes a variable type and gives ctype and rebxtype
        self.from_value =   {   'int':(c_int, val["REBX_TYPE_INT"]), 
                                'int32':(c_int, val["REBX_TYPE_INT"]), 
                                'int64':(c_int, val["REBX_TYPE_INT"]), 
                                'float':(c_double, val["REBX_TYPE_DOUBLE"]), 
                                'float64':(c_double, val["REBX_TYPE_DOUBLE"]),
                                'c_uint':(c_uint, val["REBX_TYPE_UINT32"]),
                                'c_uint8':(c_uint, val["REBX_TYPE_UINT32"]),
                                'c_uint16':(c_uint, val["REBX_TYPE_UINT32"]),
                                'c_uint32':(c_uint, val["REBX_TYPE_UINT32"]),
                                'c_uint64':(c_uint, val["REBX_TYPE_UINT32"]),
                            }

    def __getitem__(self, key):
        name = rebound.hash(key).value
        clibreboundx.rebx_get_param_node.restype = POINTER(Param)
        node = clibreboundx.rebx_get_param_node(byref(self.parent), c_char_p(key.encode('ascii')))
        if node:
            param = node.contents
            ctype, dtype = self.from_type[param.param_type]
            if param.ndim == 0:
                data = cast(param.contents, POINTER(ctype))
                return data.contents.value
            else:
                ArrayType = ctype*param.size
                data = cast(param.contents, POINTER(ArrayType))
                import numpy as np
                array = np.frombuffer(data.contents, dtype=dtype, count=param.size)
                if param.ndim > 1:
                    shape = (param.shape[i] for i in range(param.ndim))
                    array = np.reshape(array, tuple(shape))
                return array
        else:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' not found.".format(key))

    def __setitem__(self, key, value):
        if type(value) == list:
            try:
                import numpy as np
                value = np.array(value)
            except ImportError:
                raise AttributeError("Need to install numpy in order to assign lists as parameters.")
    
        if type(value).__module__ == 'numpy':   # In both cases these get a string with typename (e.g. 'float64' for numpy.float64)
            valtype = value.dtype.name          # Gets type also for ndarrays
        else:
            valtype = type(value).__name__
       
        try:
            ctype, rebxtype = self.from_value[valtype]
        except KeyError:
            raise AttributeError("REBOUNDx Error: Data type {0} for param '{1}' not supported.".format(valtype, key))
        
        clibreboundx.rebx_get_param_node.restype = POINTER(Param)
        nodeptr = clibreboundx.rebx_get_param_node(byref(self.parent), c_char_p(key.encode('ascii')))

        if nodeptr:
            param = nodeptr.contents
            if param.param_type != rebxtype:
                raise AttributeError("REBOUNDx Error: Cannot update param '{0}' with new type {1}".format(key, valtype))
            if type(value).__name__ == 'ndarray':
                if param.ndim == 0:
                    raise AttributeError("REBOUNDx Error: Cannot update scalar param '{0}' with a list or numpy array".format(key))
                else:
                    ArrayType = c_int*param.ndim
                    data = cast(param.shape, POINTER(ArrayType))
                    import numpy as np
                    shape = np.frombuffer(data.contents, dtype='int32', count=param.ndim)
                    if not np.array_equal(shape, value.shape):
                        raise AttributeError("REBOUNDx Error: Cannot update param '{0}' with list/array of different shape".format(key))
            else:
                if param.ndim > 0:
                    raise AttributeError("REBOUNDx Error: Cannot update param '{0}' (a list/array) with scalar".format(key)) 
            val = nodeptr.contents.contents
        else: # parameter not found. Make new one
            if type(value).__name__ == 'ndarray':
                import numpy as np
                clibreboundx.rebx_add_param_.restype = c_void_p
                val = clibreboundx.rebx_add_param_(byref(self.parent), c_char_p(key.encode('ascii')), rebxtype, value.ndim, value.ctypes.shape_as(c_int))
            else: # single value
                clibreboundx.rebx_add_param.restype = c_void_p
                val = clibreboundx.rebx_add_param(byref(self.parent), c_char_p(key.encode('ascii')), rebxtype)

        val = cast(val, POINTER(ctype))
        if type(value).__name__ == 'ndarray':
            if 'int' in valtype:
                value = value.astype('int32', casting='same_kind')
            import numpy as np
            memmove(val, value.ctypes.data_as(POINTER(ctype)), sizeof(ctype)*value.size) # COPIES data
        else:
            if 'int' in valtype:
                if 'uint' in valtype:
                    if value.value > 4294967295:
                        raise OverflowError("REBOUNDx Error: Unsigned integer parameters cannot exceed 2^32")
                elif value > 2147483648:
                    raise OverflowError("REBOUNDx Error: Integer parameters cannot exceed 2^31")
            val[0] = value
        
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
