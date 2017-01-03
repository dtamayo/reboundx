import rebound
from collections import MutableMapping
from .extras import Param, Effect
from . import clibreboundx
from ctypes import byref, c_double, c_int, c_int32, c_int64, c_uint, c_uint32, c_longlong, c_char_p, POINTER, cast
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
        self.ptype = {enum:ptype for (enum, ptype) in enumerate(python_types)} # get python type from enum value
        self.defaultptype = {val[rebxtype]:ptype for (rebxtype, ptype) in rebx_param_types} # get default python type from C rebx type
        self.penum = {ptype:enum for (enum, ptype) in enumerate(python_types)} # get enum from python type
        self.types = {item[0]:(item[1], item[2], item[3]) for item in type_tuples} # (dtype, ctype, rebxtype)from ptype

    def __getitem__(self, key):
        if self.verbose: print("*****")
        clibreboundx.rebx_get_param_node.restype = POINTER(Param)
        node = clibreboundx.rebx_get_param_node(byref(self.parent), c_char_p(key.encode('ascii')))
       
        if node: 
            if self.verbose: print("Parameter {0} found".format(key))
            param = node.contents
            try:
                pythontype = self.ptype[param.python_type]
            except KeyError:
                pythontype = self.defaultptype[param.param_type] # Parameter was set in C, use default

            dtype, ctype, rebxtype = self.types[pythontype]
            if self.verbose:
                print("param.python_type= {0}, pythontype = {1}".format(param.python_type, pythontype))
                print("dtype = {0}, ctype = {1}, rebxtype = {2}".format(dtype, ctype, rebxtype))

            if param.ndim == 0: 
                if self.verbose: print("Scalar")
                data = cast(param.contents, POINTER(ctype))
                try:
                    retval = data.contents.value # cases where scalar is a ctype, so we get python type back 
                    retval = pythontype(retval)
                except AttributeError:
                    if self.verbose: print("Parameter is a Structure")
                    retval = data.contents # cases where scalar is a Structure, e.g., rebound.Orbit
            else:
                if self.verbose: print("Parameter is an Array. ndim = {0}".format(param.ndim))
                if self.verbose: print("Size = {0}".format(param.size))
                ArrayType = ctype*param.size
                data = cast(param.contents, POINTER(ArrayType))
                if dtype == object: # can't use frombuffer with objects. 
                    try:
                        array = np.array([pythontype(data.contents[i]) for i in range(param.size)], dtype=dtype)
                    except TypeError: # for cases where a custom ctypes Structure can't take a class instance in __init__
                        array = np.array([data.contents[i] for i in range(param.size)], dtype=dtype)
                    #Makes copies of the ctypes.Structure instances, but fields share memory with C.
                else:
                    array = np.frombuffer(data.contents, dtype=dtype, count=param.size)
                if param.ndim > 1:
                    shape = [param.shape[i] for i in range(param.ndim)]
                    array = np.reshape(array, shape)
                    if self.verbose: print("Array shape = {0}".format(shape))
                retval = array
        else:
            raise AttributeError("REBOUNDx Error: Parameter '{0}' not found.".format(key))
        
        if self.verbose: print("*****")
        return retval

    def __setitem__(self, key, value):
        if type(value) == list:
            raise AttributeError("Can't assign lists as params. Use numpy arrays instead. See documentation.")
            '''
            value = np.array(value, dtype=object) 
            ''' 
        tempval = value
        try:
            while True: # keep going through nested arrays until we get to a real entry (that we can't index)
                tempval = tempval[0]
                if type(tempval).__name__ == 'void': # this is the type you get if you don't use dtype=object
                                                     # (i.e. when np autogenerates a dtype (which we can't detect)
                    raise AttributeError('When setting a numpy array of objects, have to set dtype=object, '
                            'e.g., ps[1].params["array"] = numpy.array([obj1, obj2, obj3], dtype=object)')
        except:
            pythontype = type(tempval)
            if self.verbose: print("1st array item (to check type) = {0}, type = {1}".format(tempval, pythontype))

        try:
            dtype, ctype, rebxtype = self.types[pythontype]
            if self.verbose: print("ctype = {0}, rebxtype = {1}".format(ctype, rebxtype))
        except KeyError:
            raise AttributeError("REBOUNDx Error: Data type {0} for param '{1}' "
            "not supported.".format(pythontype, key))
        
        clibreboundx.rebx_get_param_node.restype = POINTER(Param)
        nodeptr = clibreboundx.rebx_get_param_node(byref(self.parent), c_char_p(key.encode('ascii')))

        if nodeptr:
            if self.verbose: print("Parameter {0} found.".format(key))
            param = nodeptr.contents
            if param.param_type != rebxtype:
                raise AttributeError("REBOUNDx Error: Can't update param '{0}' with "
                "new type {1}".format(key, pythontype))
            if type(value).__name__ == 'ndarray':
                if self.verbose: print("Param is an array")
                if param.ndim == 0:
                    raise AttributeError("REBOUNDx Error: Cannot update scalar param '{0}' with "
                    "a list or numpy array".format(key))
                else:
                    ArrayType = c_int*param.ndim # shape is an integer array, so we recast it to ctypes array
                    data = cast(param.shape, POINTER(ArrayType))
                    shape = np.frombuffer(data.contents, dtype='int32', count=param.ndim)
                    if self.verbose: print("Array shape = {0}".format(shape))
                    if not np.array_equal(shape, value.shape):
                        raise AttributeError("REBOUNDx Error: Cannot update param '{0}' with a "
                        "list/array of different shape".format(key))
            else:
                if param.ndim > 0:
                    raise AttributeError("REBOUNDx Error: Cannot update param '{0}' "
                    "(a list/array) with scalar".format(key)) 
            val = nodeptr.contents.contents
        else: # parameter not found. Make new one
            if self.verbose: print("Parameter {0} not found. Making new one".format(key))
            if type(value) == np.ndarray:
                if self.verbose: print("Param is an array of dimension {0}".format(value.ndim))
                clibreboundx.rebx_add_param_node.restype = POINTER(Param)
                nodeptr = clibreboundx.rebx_add_param_node(
                        byref(self.parent), c_char_p(key.encode('ascii')), rebxtype, 
                        value.ndim, value.ctypes.shape_as(c_int))
            else: # single value
                if self.verbose: print("Param is a scalar")
                clibreboundx.rebx_add_param_node.restype = POINTER(Param)
                nodeptr = clibreboundx.rebx_add_param_node(
                        byref(self.parent), c_char_p(key.encode('ascii')), rebxtype, 
                        0, None)
            nodeptr.contents.python_type = c_int(self.penum[pythontype])
            if self.verbose:
                print("param.python_type = {0}, pythontype = {1}".format(nodeptr.contents.python_type, pythontype))

        val = cast(nodeptr.contents.contents, POINTER(ctype))
        if type(value) == np.ndarray:
            ArrayType = ctype*value.size
            ctypesarray = ArrayType(*(i for i in value.flatten())) # need to flatten to 1D (final format anyway)
                                                                   # so we don't populate ctypesarray with 
                                                                   # arrays for multidim arrays
            memmove(val, ctypesarray, sizeof(ctype)*value.size)    # COPIES data so that C owns the memory
        else:
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
