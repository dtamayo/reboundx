import rebound
from collections import MutableMapping
from .extras import Param, Effect
from . import clibreboundx
from ctypes import byref, c_double, c_int, c_int32, c_int64, c_uint, c_uint32, c_char_p, POINTER, cast
from ctypes import c_void_p, memmove, sizeof, addressof
import numpy as np

class Params(MutableMapping):
    def __init__(self, parent):
        self.verbose = 0 # set to 1 to diagnose problems
        self.parent = parent

        # Check rebx instance is attached, otherwise memory isn't allocated.
        if parent.__class__ == rebound.particle.Particle: 
            if not parent._sim.contents.extras:
                raise AttributeError("Need to attach reboundx.Extras instance to simulation before setting"
                        " particle params.")

        #### Update here to add new types #####
        
        rebx_param_types = ["REBX_TYPE_DOUBLE", 
                            "REBX_TYPE_INT", 
                            "REBX_TYPE_UINT32", 
                            "REBX_TYPE_ORBIT",
                           ] # In same order as C enum
        val = {param_type:value for (value, param_type) in enumerate(rebx_param_types)}

        # takes rebxtype and gives ctype and dtype for numpy arrays
        self.from_type =    {val["REBX_TYPE_INT"]:(c_int, 'int32'), # REBOUNDx uses 32 bit int for all int params 
                             val["REBX_TYPE_DOUBLE"]:(c_double, 'float'),
                             val["REBX_TYPE_UINT32"]:(c_uint32, 'uint32'),
                             val["REBX_TYPE_ORBIT"]:(rebound.Orbit, 'object'),
                            }

        # takes a variable type and gives ctype and rebxtype
        self.from_value =   {'int':(c_int, val["REBX_TYPE_INT"]), 
                             'int32':(c_int, val["REBX_TYPE_INT"]), 
                             'int64':(c_int, val["REBX_TYPE_INT"]), 
                             'float':(c_double, val["REBX_TYPE_DOUBLE"]), 
                             'float64':(c_double, val["REBX_TYPE_DOUBLE"]),
                             'c_uint':(c_uint, val["REBX_TYPE_UINT32"]),
                             'c_uint8':(c_uint, val["REBX_TYPE_UINT32"]),
                             'c_uint16':(c_uint, val["REBX_TYPE_UINT32"]),
                             'c_uint32':(c_uint, val["REBX_TYPE_UINT32"]),
                             'c_uint64':(c_uint, val["REBX_TYPE_UINT32"]),
                             'Orbit':(rebound.Orbit, val["REBX_TYPE_ORBIT"]),
                            }

    def __getitem__(self, key):
        if self.verbose: print("*****")

        clibreboundx.rebx_get_param_node.restype = POINTER(Param)
        node = clibreboundx.rebx_get_param_node(byref(self.parent), c_char_p(key.encode('ascii')))
        
        if node: 
            if self.verbose: print("Parameter {0} found".format(key))
            param = node.contents
            ctype, dtype = self.from_type[param.param_type]
            if self.verbose:
                print("ctype = {0}, param_type = {2}".format(ctype, dtype, param.param_type))

            if param.ndim == 0: 
                if self.verbose: print("Scalar")
                data = cast(param.contents, POINTER(ctype))
                try:
                    retval = data.contents.value # cases where scalar is a ctype, so we get python type back 
                except AttributeError:
                    if self.verbose: print("Parameter is a Structure")
                    retval = data.contents # cases where scalar is a Structure, e.g., rebound.Orbit
            else:
                if self.verbose: print("Parameter is an Array. ndim = {0}".format(param.ndim))
                try:
                    import numpy as np
                except ImportError:
                    raise AttributeError("Need to install numpy in order to use array parameters.")
                
                if self.verbose: print("Size = {0}".format(param.size))
                ArrayType = ctype*param.size
                data = cast(param.contents, POINTER(ArrayType))

                if self.verbose: print("dtype = {0}".format(dtype))
                if dtype == 'object': # can't use frombuffer with objects. 
                    array = np.array([data.contents[i] for i in range(param.size)], dtype=dtype)
                    #Makes copies of the ctypes.Structure instances, but fields share memory with C.
                else:
                    array = np.frombuffer(data.contents, dtype=dtype, count=param.size)# shares memory with C

                if param.ndim > 1:
                    shape = (param.shape[i] for i in range(param.ndim))
                    shape = tuple(shape)
                    array = np.reshape(array, shape)
                    print("Array shape = {0}".format(shape))
                retval = array
        else:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' not found.".format(key))
        
        if self.verbose: print("*****")
        return retval

    def __setitem__(self, key, value):
        if type(value) == list:
            try:
                import numpy as np
                value = np.array(value, dtype=object) 
                # with list of class instances, need dtype=object so we can detect type later.  
                # Otherwise numpy autogenerates dtype from all fields, which we can't use.
            except ImportError:
                raise AttributeError("Need to install numpy in order to assign lists as parameters.")
     
        tempval = value
        try:
            while True: # keep going through nested arrays until we get to a real entry (that we can't index)
                tempval = tempval[0]
                if type(tempval).__name__ == 'void': # this is the type you get if you don't use dtype=object
                                                     # (i.e. when np autogenerates a dtype (which we can't detect)
                    raise AttributeError('When setting a numpy array of objects, have to set dtype=object, '
                            'e.g., ps[1].params["array"] = numpy.array([obj1, obj2, obj3], dtype=object)')
        except:
            valtype = type(tempval).__name__

        try:
            ctype, rebxtype = self.from_value[valtype]
        except KeyError:
            raise AttributeError("REBOUNDx Error: Data type {0} for param '{1}' "
            "not supported.".format(valtype, key))
        
        clibreboundx.rebx_get_param_node.restype = POINTER(Param)
        nodeptr = clibreboundx.rebx_get_param_node(byref(self.parent), c_char_p(key.encode('ascii')))

        if nodeptr:
            param = nodeptr.contents
            if param.param_type != rebxtype:
                raise AttributeError("REBOUNDx Error: Can't update param '{0}' with "
                "new type {1}".format(key, valtype))
            if type(value).__name__ == 'ndarray':
                if param.ndim == 0:
                    raise AttributeError("REBOUNDx Error: Cannot update scalar param '{0}' with "
                    "a list or numpy array".format(key))
                else:
                    ArrayType = c_int*param.ndim # shape is an integer array, so we recast it to ctypes array
                    data = cast(param.shape, POINTER(ArrayType))
                    import numpy as np
                    shape = np.frombuffer(data.contents, dtype='int32', count=param.ndim)
                    if not np.array_equal(shape, value.shape):
                        raise AttributeError("REBOUNDx Error: Cannot update param '{0}' with a "
                        "list/array of different shape".format(key))
            else:
                if param.ndim > 0:
                    raise AttributeError("REBOUNDx Error: Cannot update param '{0}' "
                    "(a list/array) with scalar".format(key)) 
            val = nodeptr.contents.contents
        else: # parameter not found. Make new one
            if type(value).__name__ == 'ndarray':
                import numpy as np
                clibreboundx.rebx_add_param_.restype = c_void_p
                val = clibreboundx.rebx_add_param_(
                        byref(self.parent), c_char_p(key.encode('ascii')), rebxtype, 
                        value.ndim, value.ctypes.shape_as(c_int))
            else: # single value
                clibreboundx.rebx_add_param.restype = c_void_p
                val = clibreboundx.rebx_add_param(byref(self.parent), c_char_p(key.encode('ascii')), rebxtype)

        val = cast(val, POINTER(ctype))
        if type(value).__name__ == 'ndarray':
            ArrayType = ctype*value.size
            ctypesarray = ArrayType(*(i for i in value.flatten())) # need to flatten to 1D (final format anyway)
                                                                   # so we don't populate ctypesarray with 
                                                                   # arrays for multidim arrays
            memmove(val, ctypesarray, sizeof(ctype)*value.size)    # COPIES data so that C owns the memory
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
