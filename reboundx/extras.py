from . import clibreboundx
from ctypes import Structure, c_double, POINTER, c_int, c_uint, c_long, c_ulong, c_void_p, c_char_p, CFUNCTYPE, byref, c_uint32, c_uint, cast, c_char
import rebound
import reboundx
import warnings

INTEGRATORS = {"implicit_midpoint": 0, "rk4":1, "euler": 2, "rk2": 3, "none": -1}

REBX_TIMING = {"pre":-1, "post":1}
REBX_FORCE_TYPE = {"none":0, "pos":1, "vel":2}
REBX_OPERATOR_TYPE = {"none":0, "updater":1, "recorder":2}

REBX_BINARY_WARNINGS = [
        ("REBOUNDx Error: Cannot read binary file. Check filename and file contents.", 1),
        ("REBOUNDx Error: Binary file was corrupt. Could not read.", 2),
        ("REBOUNDx Warning: Binary file was saved with a different version of REBOUNDx. Binary format might have changed.", 4),
        ("REBOUNDx Warning: At least one parameter in the binary file was not loaded. Check simulation.", 8),
        ("REBOUNDx Warning: At least one particle's parameters in the binary file were not loaded. Check simulation.", 16),
        ("REBOUNDx Warning: At least one effect and its parameters were not loaded from the binary file. Check simulation.", 32),
        ("REBOUNDx Warning: At least one field in the binary field was not recognized, and not loaded. Probably binary was created with more recent REBOUNDx version than you are using.", 64)]

class Extras(Structure):
    """
    Main object used for all REBOUNDx operations, tied to a particular REBOUND simulation.
    This is an abstraction of the C struct rebx_extras, with all the C convenience functions
    and functions for adding effects implemented as methods of the class.  
    The fastest way to understand it is to follow the examples at :ref:`ipython_examples`.  
    """
    
    def __new__(cls, sim, filename=None):
        # Create simulation
        if filename==None:
            # Create a new rebx instance
            rebx = super(Extras,cls).__new__(cls)
            clibreboundx.rebx_initialize(byref(sim), byref(rebx))
            clibreboundx.rebx_register_default_params(byref(self))
        else:
            pass
            '''
            # Recreate exisitng simulation 
            sa = SimulationArchive(filename,process_warnings=False)
            sim = super(Simulation,cls).__new__(cls)
            clibrebound.reb_init_simulation(byref(sim))
            w = sa.warnings # warnings will be appended to previous warnings (as to not repeat them) 
            clibrebound.reb_create_simulation_from_simulationarchive_with_messages(byref(sim),byref(sa),c_int(snapshot),byref(w))
            for majorerror, value, message in BINARY_WARNINGS:
                if w.value & value:
                    if majorerror:
                        raise RuntimeError(message)
                    else:  
                        # Just a warning
                        warnings.warn(message, RuntimeWarning)
            '''
        if not hasattr(sim, "_extras_ref"): # if REBOUNDx wasn't already attached, check for warnings in case additional_forces or ptm were already set.
            sim.process_messages()
        sim._extras_ref = rebx # add a reference to this instance in sim to make sure it's not garbage collected
        rebx.custom_effects = {} # dictionary to keep references to custom effects so they don't get garbage collected
        return rebx 

    def __init__(self, sim, filename=None):
        pass

    @classmethod
    def from_file(cls, sim, filename):
        """
        Loads REBOUNDx effects along with effect and particle parameters from a binary file.
        
        Arguments
        ---------
        filename : str
            Filename of the binary file.
        
        Returns
        ------- 
        A reboundx.Extras object.
        
        """
        w = c_int(0)
        clibreboundx.rebx_init.restype = POINTER(Extras)
        extrasp = clibreboundx.rebx_init(byref(sim))
        clibreboundx.rebx_create_extras_from_binary_with_messages(extrasp, c_char_p(filename.encode("ascii")),byref(w))
        if (extrasp is None) or (w.value & 1):     # Major error
            raise ValueError(REBX_BINARY_WARNINGS[0])
        for message, value in REBX_BINARY_WARNINGS:  # Just warnings
            if w.value & value and value!=1:
                warnings.warn(message, RuntimeWarning)
        extras = extrasp.contents
        sim.save_messages = 1 # Warnings will be checked within python
        return extras 

    def __del__(self):
        if self._b_needsfree_ == 1:
            clibreboundx.rebx_free_pointers(byref(self))

    def detach(self):
        clibreboundx.rebx_detach(self._sim)

    @property
    def integrator(self):
        """
        Get or set the intergrator module.

        Available integrators are:

        - ``'implicit_midpoint'`` (default)
        
        Check the online documentation for a full description of each of the integrators. 
        """
        i = self._integrator
        for name, _i in INTEGRATORS.items():
            if i==_i:
                return name
        return i
    @integrator.setter
    def integrator(self, value):
        if isinstance(value, int):
            self._integrator = c_int(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in INTEGRATORS: 
                self._integrator = INTEGRATORS[value]
            else:
                raise ValueError("Warning. Integrator not found.")
    
    #######################################
    # Functions for manipulating REBOUNDx effects
    #######################################

    def register_param(self, name, param_type):
        type_enum = REBX_C_PARAM_TYPES[param_type] 
        clibreboundx.rebx_register_param(byref(self), c_char_p(name.encode('ascii')), c_int(type_enum))
        self._sim.contents.process_messages()

    def load_force(self, name):
        clibreboundx.rebx_load_force.restype = POINTER(Force)
        ptr = clibreboundx.rebx_load_force(byref(self), c_char_p(name.encode('ascii')))
        self._sim.contents.process_messages()
        return ptr.contents
    
    def create_force(self, name):
        clibreboundx.rebx_create_force.restype = POINTER(Force)
        ptr = clibreboundx.rebx_create_force(byref(self), c_char_p(name.encode('ascii')))
        self._sim.contents.process_messages()
        return ptr.contents

    def load_operator(self, name):
        clibreboundx.rebx_load_operator.restype = POINTER(Operator)
        ptr = clibreboundx.rebx_load_operator(byref(self), c_char_p(name.encode('ascii')))
        self._sim.contents.process_messages()
        return ptr.contents

    def create_operator(self, name):
        clibreboundx.rebx_create_operator.restype = POINTER(Operator)
        ptr = clibreboundx.rebx_create_operator(byref(self), c_char_p(name.encode('ascii')))
        self._sim.contents.process_messages()
        return ptr.contents

    def add_force(self, force):
        clibreboundx.rebx_add_force(byref(self), byref(force))
        self._sim.contents.process_messages()

    def add_operator(self, operator, dt_fraction=None, timing="post", name=""):
        if dt_fraction is None:
            clibreboundx.rebx_add_operator(byref(self), byref(operator))
        else:
            timingint = REBX_TIMING[timing]
            clibreboundx.rebx_add_operator_step(byref(self), byref(operator), c_double(dt_fraction), c_int(timingint), c_char_p(name.encode('ascii')))
        self._sim.contents.process_messages()

    def get_force(self, name):
        clibreboundx.rebx_get_force.restype = POINTER(Force)
        ptr = clibreboundx.rebx_get_force(byref(self), c_char_p(name.encode('ascii')))
        if ptr:
            return ptr.contents
        else:
            raise AttributeError("REBOUNDx Error: Force {0} passed to rebx.get_force not found.".format(name))

    def get_operator(self, name):
        clibreboundx.rebx_get_operator.restype = POINTER(Operator)
        ptr = clibreboundx.rebx_get_operator(byref(self), c_char_p(name.encode('ascii')))
        self._sim.contents.process_messages()
        if ptr:
            return ptr.contents
        else:
            raise AttributeError("REBOUNDx Error: Force {0} passed to rebx.get_operator not found.".format(name))

    def remove_force(self, force):
        if not isinstance(force, reboundx.extras.Force):
            raise TypeError("REBOUNDx Error: Object passed to rebx.remove_force is not a reboundx.Force instance.")
        success = clibreboundx.rebx_remove_force(byref(self), byref(force))
        if not success:
            raise AttributeError("REBOUNDx Error: Force {0} passed to rebx.remove_force not found in simulation.")

    def remove_operator(self, operator):
        if not isinstance(operator, reboundx.extras.Operator):
            raise TypeError("REBOUNDx Error: Object passed to rebx.remove_operator is not a reboundx.Operator instance.")
        success = clibreboundx.rebx_remove_operator(byref(self), byref(operator))
        if not success:
            raise AttributeError("REBOUNDx Error: Operator {0} passed to rebx.remove_operator not found in simulation.")

    #######################################
    # Input/Output Routines
    #######################################
    def save(self, filename):
        """
        Save the entire REBOUND simulation to a binary file.
        """
        clibreboundx.rebx_output_binary(byref(self), c_char_p(filename.encode("ascii")))

    #######################################
    # Convenience Functions
    #######################################

    def rad_calc_beta(self, G, c, source_mass, source_luminosity, radius, density, Q_pr):
        clibreboundx.rebx_rad_calc_beta.restype = c_double
        return clibreboundx.rebx_rad_calc_beta(c_double(G), c_double(c), c_double(source_mass), c_double(source_luminosity), c_double(radius), c_double(density), c_double(Q_pr))

    def rad_calc_particle_radius(self, G, c, source_mass, source_luminosity, beta, density, Q_pr):
        clibreboundx.rebx_rad_calc_particle_radius.restype = c_double
        return clibreboundx.rebx_rad_calc_particle_radius(c_double(G), c_double(c), c_double(source_mass), c_double(source_luminosity), c_double(beta), c_double(density), c_double(Q_pr))

    def central_force_Acentral(self, p, primary, pomegadot, gamma):
        clibreboundx.rebx_central_force_Acentral.restype = c_double
        Acentral = clibreboundx.rebx_central_force_Acentral(p, primary, c_double(pomegadot), c_double(gamma))
        self._sim.contents.process_messages()
        return Acentral

    # Hamiltonian calculation functions
    def gr_hamiltonian(self, sim, params):
        clibreboundx.rebx_gr_hamiltonian.restype = c_double
        return clibreboundx.rebx_gr_hamiltonian(byref(sim), byref(params))
    
    def gr_potential_hamiltonian(self, sim, params):
        clibreboundx.rebx_gr_potential_hamiltonian.restype = c_double
        return clibreboundx.rebx_gr_potential_hamiltonian(byref(sim), byref(params))
    
    def gr_full_hamiltonian(self, sim, params):
        clibreboundx.rebx_gr_full_hamiltonian.restype = c_double
        return clibreboundx.rebx_gr_full_hamiltonian(byref(sim), byref(params))
    
    def tides_precession_hamiltonian(self, sim, params):
        clibreboundx.rebx_tides_precession_hamiltonian.restype = c_double
        return clibreboundx.rebx_tides_precession_hamiltonian(byref(sim), byref(params))

    def central_force_hamiltonian(self, sim):
        clibreboundx.rebx_central_force_hamiltonian.restype = c_double
        return clibreboundx.rebx_central_force_hamiltonian(byref(sim))
    
    def gravitational_harmonics_hamiltonian(self, sim):
        clibreboundx.rebx_gravitational_harmonics_hamiltonian.restype = c_double
        return clibreboundx.rebx_gravitational_harmonics_hamiltonian(byref(sim))
    
#################################################
# Generic REBOUNDx definitions
#################################################


class Param(Structure): # need to define fields afterward because of circular ref in linked list
    pass    
Param._fields_ =  [ ("name", c_char_p),
                    ("type", c_int),
                    ("value", c_void_p)]

class Node(Structure): # need to define fields afterward because of circular ref in linked list
    pass    
Node._fields_ =  [  ("object", c_void_p),
                    ("next", POINTER(Node))]

class Operator(Structure):
    @property
    def operator_type(self):
        return self._operator_type

    @operator_type.setter
    def operator_type(self, value):
        self._operator_type = REBX_OPERATOR_TYPE[value.lower()]

    @property
    def step(self):
        return self._step

    @step.setter
    def step(self, func):
        self._sfp = STEPFUNCPTR(func) # keep a reference to func so it doesn't get garbage collected
        self._step = self._sfp

    @property 
    def params(self):
        params = Params(self)
        return params

STEPFUNCPTR = CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(Operator), c_double)

Operator._fields_ = [   ("name", c_char_p),
                        ("ap", POINTER(Node)),
                        ("_sim", POINTER(rebound.Simulation)),
                        ("_operator_type", c_int),
                        ("_step", STEPFUNCPTR)]
class Force(Structure):
    @property
    def force_type(self):
        return self._force_type

    @force_type.setter
    def force_type(self, value):
        self._force_type = REBX_FORCE_TYPE[value.lower()]

    @property
    def update_accelerations(self):
        return self._update_accelerations

    @update_accelerations.setter
    def update_accelerations(self, func):
        self._ffp = FORCEFUNCPTR(func) # keep a reference to func so it doesn't get garbage collected
        self._update_accelerations = self._ffp

    @property 
    def params(self):
        params = Params(self)
        return params

FORCEFUNCPTR = CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(Force), POINTER(rebound.Particle), c_int)

Force._fields_ = [  ("name", c_char_p),
                    ("ap", POINTER(Node)),
                    ("_sim", POINTER(rebound.Simulation)),
                    ("_force_type", c_int),
                    ("_update_accelerations", FORCEFUNCPTR)]

# Need to put fields after class definition because of self-referencing
Extras._fields_ =  [("_sim", POINTER(rebound.Simulation)),
                    ("_additional_forces", POINTER(Node)),
                    ("_pre_timestep_modifications", POINTER(Node)),
                    ("_post_timestep_modifications", POINTER(Node)),
                    ("_registered_params", POINTER(Node)),
                    ("_allocated_forces", POINTER(Node)),
                    ("_allocated_operators", POINTER(Node)),
                    ("_integrator", c_int)]

# This list keeps pairing from C rebx_param_type enum to ctypes type 1-to-1. Derive the required mappings from it
REBX_C_TO_CTYPES = [["REBX_TYPE_NONE", None], ["REBX_TYPE_DOUBLE", c_double], ["REBX_TYPE_INT",c_int], ["REBX_TYPE_POINTER", c_void_p], ["REBX_TYPE_FORCE", Force]]
REBX_CTYPES = {} # maps int value of rebx_param_type enum to ctypes type
REBX_C_PARAM_TYPES = {} # maps string of rebx_param_type enum to int
for i, pair in enumerate(REBX_C_TO_CTYPES):
    REBX_CTYPES[i] = pair[1]
    REBX_C_PARAM_TYPES[pair[0]] = i

from .params import Params
