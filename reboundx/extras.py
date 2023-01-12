from . import clibreboundx
from ctypes import Structure, c_double, POINTER, c_int, c_uint, c_long, c_ulong, c_void_p, c_char_p, CFUNCTYPE, byref, c_uint32, c_uint, cast, c_char, pointer
import rebound
import reboundx
import warnings

integrators = {"implicit_midpoint": 0, "rk4":1, "euler": 2, "rk2": 3, "none": -1}

REBX_TIMING = {"pre":-1, "post":1}
REBX_FORCE_TYPE = {"none":0, "pos":1, "vel":2}
REBX_OPERATOR_TYPE = {"none":0, "updater":1, "recorder":2}

REBX_BINARY_WARNINGS = [
    (True, 1, "REBOUNDx: Cannot open binary file. Check filename."),
    (True, 2, "REBOUNDx: Binary file is unreadable. Please open an issue on Github mentioning the version of REBOUND and REBOUNDx you are using and include the binary file."),
    (True, 4, "REBOUNDx: Ran out of system memory."),
    (True, 8, "REBOUNDx: REBOUNDx structure couldn't be loaded."),
    (True, 16, "REBOUNDx: At least one registered parameter was not loaded. This typically indicates the binary is corrupt or was saved with an incompatible version to the current one being used."),
    (False, 32, "REBOUNDx: At least one force or operator parameter was not loaded from the binary file. This typically indicates the binary is corrupt or was saved with an incompatible version to the current one being used."),
    (False, 64, "REBOUNDx: At least one particle's parameters were not loaded from the binary file."),
    (False, 128, "REBOUNDx: At least one force was not loaded from the binary file. If binary was created with a newer version of REBOUNDx, a particular force may not be implemented in your current version of REBOUNDx."),
    (False, 256, "REBOUNDx: At least one operator was not loaded from the binary file. If binary was created with a newer version of REBOUNDx, a particular force may not be implemented in your current version of REBOUNDx."),
    (False, 512, "REBOUNDx: At least one operator step was not loaded from the binary file."),
    (False, 1024, "REBOUNDx: At least one force was not added to the simulation. If binary was created with a newer version of REBOUNDx, a particular force may not be implemented in your current version of REBOUNDx."),
    (False, 2048, "REBOUNDx: Unknown field found in binary file. Any unknown fields not loaded.  This can happen if the binary was created with a later version of REBOUNDx than the one used to read it."),
    (False, 4096, "REBOUNDx: Unknown list in the REBOUNDx structure wasn't loaded. This can happen if the binary was created with a later version of REBOUNDx than the one used to read it."),
    (False, 8192, "REBOUNDx: The value of at least one parameter was not loaded. This can happen if a custom structure was added by the user as a parameter. See Parameters.ipynb jupyter notebook example."),
    (False,16384, "REBOUNDx: Binary file was saved with a different version of REBOUNDx. Binary format might have changed. Check that effects and parameters are loaded as expected.")
]

class Extras(Structure):
    """
    Main object used for all REBOUNDx operations, tied to a particular REBOUND simulation.
    This is an abstraction of the C struct rebx_extras, with all the C convenience functions
    and functions for adding effects implemented as methods of the class.
    The fastest way to understand it is to follow the examples at :ref:`ipython_examples`.
    """

    def __new__(cls, sim, filename=None):
        rebx = super(Extras,cls).__new__(cls)
        return rebx

    def __init__(self, sim, filename=None):
        sim._extras_ref = self # add a reference to this instance in sim to make sure it's not garbage collected_
        clibreboundx.rebx_initialize(byref(sim), byref(self))
        # Create simulation
        if filename==None:
            # Create a new rebx instance
           clibreboundx.rebx_register_default_params(byref(self))
        else:
            # Recreate existing simulation.
            # Load registered parameters from binary
            w = c_int(0)
            clibreboundx.rebx_init_extras_from_binary(byref(self), c_char_p(filename.encode('ascii')), byref(w))
            for majorerror, value, message in REBX_BINARY_WARNINGS:
                if w.value & value:
                    if majorerror:
                        raise RuntimeError(message)
                    else:
                        warnings.warn(message, RuntimeWarning)
        self.process_messages()

    def __del__(self):
        if self._b_needsfree_ == 1:
            clibreboundx.rebx_free_pointers(byref(self))

    def detach(self, sim):
        sim._extras_ref = None # remove reference to rebx so it can be garbage collected
        clibreboundx.rebx_detach(byref(sim), byref(self))

    #######################################
    # Functions for manipulating REBOUNDx effects
    #######################################

    def register_param(self, name, param_type):
        type_enum = REBX_C_PARAM_TYPES[param_type]
        clibreboundx.rebx_register_param(byref(self), c_char_p(name.encode('ascii')), c_int(type_enum))
        self.process_messages()

    def load_force(self, name):
        clibreboundx.rebx_load_force.restype = POINTER(Force)
        ptr = clibreboundx.rebx_load_force(byref(self), c_char_p(name.encode('ascii')))
        self.process_messages()
        return ptr.contents

    def create_force(self, name):
        clibreboundx.rebx_create_force.restype = POINTER(Force)
        ptr = clibreboundx.rebx_create_force(byref(self), c_char_p(name.encode('ascii')))
        self.process_messages()
        return ptr.contents

    def load_operator(self, name):
        clibreboundx.rebx_load_operator.restype = POINTER(Operator)
        ptr = clibreboundx.rebx_load_operator(byref(self), c_char_p(name.encode('ascii')))
        self.process_messages()
        return ptr.contents

    def create_operator(self, name):
        clibreboundx.rebx_create_operator.restype = POINTER(Operator)
        ptr = clibreboundx.rebx_create_operator(byref(self), c_char_p(name.encode('ascii')))
        self.process_messages()
        return ptr.contents

    def add_force(self, force):
        if not isinstance(force, reboundx.extras.Force):
            raise TypeError("REBOUNDx Error: Object passed to rebx.add_force is not a reboundx.Force instance.")
        clibreboundx.rebx_add_force(byref(self), byref(force))
        self.process_messages()

    def add_operator(self, operator, dtfraction=None, timing="post"):
        if not isinstance(operator, reboundx.extras.Operator):
            raise TypeError("REBOUNDx Error: Object passed to rebx.add_operator is not a reboundx.Operator instance.")
        if dtfraction is None:
            clibreboundx.rebx_add_operator(byref(self), byref(operator))
        else:
            timingint = REBX_TIMING[timing]
            clibreboundx.rebx_add_operator_step(byref(self), byref(operator), c_double(dtfraction), c_int(timingint))
        self.process_messages()

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
        self.process_messages()
        if ptr:
            return ptr.contents
        else:
            raise AttributeError("REBOUNDx Error: Operator {0} passed to rebx.get_operator not found.".format(name))

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
        self.process_messages()

    #######################################
    # Effect Specific Functions
    #######################################
    
    def initialize_spin_ode(self, force):
        clibreboundx.rebx_spin_initialize_ode.restype = None
        return clibreboundx.rebx_spin_initialize_ode(byref(self), byref(force))

    def rad_calc_beta(self, G, c, source_mass, source_luminosity, radius, density, Q_pr):
        clibreboundx.rebx_rad_calc_beta.restype = c_double
        return clibreboundx.rebx_rad_calc_beta(c_double(G), c_double(c), c_double(source_mass), c_double(source_luminosity), c_double(radius), c_double(density), c_double(Q_pr))

    def rad_calc_particle_radius(self, G, c, source_mass, source_luminosity, beta, density, Q_pr):
        clibreboundx.rebx_rad_calc_particle_radius.restype = c_double
        return clibreboundx.rebx_rad_calc_particle_radius(c_double(G), c_double(c), c_double(source_mass), c_double(source_luminosity), c_double(beta), c_double(density), c_double(Q_pr))

    def central_force_Acentral(self, p, primary, pomegadot, gamma):
        clibreboundx.rebx_central_force_Acentral.restype = c_double
        Acentral = clibreboundx.rebx_central_force_Acentral(p, primary, c_double(pomegadot), c_double(gamma))
        self.process_messages()
        return Acentral

    # Hamiltonian calculation functions
    def gr_full_hamiltonian(self, force):
        clibreboundx.rebx_gr_full_hamiltonian.restype = c_double
        return clibreboundx.rebx_gr_full_hamiltonian(byref(self), byref(force))

    def gr_hamiltonian(self, force):
        clibreboundx.rebx_gr_hamiltonian.restype = c_double
        return clibreboundx.rebx_gr_hamiltonian(byref(self), byref(force))

    # Potential calculation functions
    def gr_potential_potential(self, force):
        clibreboundx.rebx_gr_potential_potential.restype = c_double
        return clibreboundx.rebx_gr_potential_potential(byref(self), byref(force))

    def tides_constant_time_lag_potential(self, force):
        clibreboundx.rebx_tides_constant_time_lag_potential.restype = c_double
        return clibreboundx.rebx_tides_constant_time_lag_potential(byref(self), byref(force))

    def tides_spin_energy(self):
        clibreboundx.rebx_tides_spin_energy.restype = c_double
        return clibreboundx.rebx_tides_spin_energy(byref(self))

    def central_force_potential(self):
        clibreboundx.rebx_central_force_potential.restype = c_double
        return clibreboundx.rebx_central_force_potential(byref(self))

    def gravitational_harmonics_potential(self):
        clibreboundx.rebx_gravitational_harmonics_potential.restype = c_double
        return clibreboundx.rebx_gravitational_harmonics_potential(byref(self))

    # Functions to help with rotations

    def rotate_simulation(self, q):
        """
        Rotates the simulation with the passed rebound.Rotation object. Analogous to Simulation.rotate, except one should use
        this function to not only rotate the orbits according to q, but also all the spatial REBOUNDx parameters (e.g., spins). 
        See SpinsIntro.ipynb
        """
        if not isinstance(q, rebound.Rotation):
            raise NotImplementedError
        clibreboundx.rebx_simulation_irotate.restype = None
        clibreboundx.rebx_simulation_irotate(byref(self), q)
  
    def spin_angular_momentum(self):
        """
        Returns a list of the three (x,y,z) components of the spin angular momentum of all particles in the simulation with
        moment of inertia (I) and spin angular frequency vector (Omega) parameters set.
        """
        clibreboundx.rebx_tools_spin_angular_momentum.restype = rebound.Vec3d
        L = clibreboundx.rebx_tools_spin_angular_momentum(byref(self))
        return [L.x, L.y, L.z]

    def process_messages(self):
        try:
            self._sim.contents.process_messages()
        except ValueError: # _sim is NULL
            raise AttributeError("REBOUNDx Error: The Simulation instance REBOUNDx was attached to no longer exists. This can happen if the Simulation instance goes out of scope or otherwise gets garbage collected.")
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
    def step_function(self):
        return self._step_function

    @step_function.setter
    def step_function(self, func):
        self._sfp = STEPFUNCPTR(func) # keep a reference to func so it doesn't get garbage collected
        self._step_function = self._sfp

    def step(self, sim, dt):
        self._step_function(byref(sim), byref(self), dt)

    @property
    def params(self):
        params = Params(self)
        return params

STEPFUNCPTR = CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(Operator), c_double)

Operator._fields_ = [   ("name", c_char_p),
                        ("ap", POINTER(Node)),
                        ("_sim", POINTER(rebound.Simulation)),
                        ("_operator_type", c_int),
                        ("_step_function", STEPFUNCPTR)]
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
                    ("_allocated_operators", POINTER(Node))]

class Interpolator(Structure):
    def __new__(cls, rebx, times, values, interpolation):
        interp = super(Interpolator, cls).__new__(cls)
        return interp

    def __init__(self, rebx, times, values, interpolation="spline"):
        try:
            Nvalues = len(times)
            Nvalues2 = len(values)
        except:
            raise TypeError("REBOUNDx Error: Times and values passed to Interpolator must be lists or arrays")
        if Nvalues != Nvalues2:
            raise ValueError("REBOUNDx Error: Times and values must be same length)")

        interpolation = interpolation.lower()
        if interpolation in INTERPOLATION_TYPE:
            interp = INTERPOLATION_TYPE[interpolation]
        else:
            raise ValueError("REBOUNDx Error: Interpolation type not supported")

        DblArr = c_double * Nvalues
        clibreboundx.rebx_init_interpolator(byref(rebx), byref(self), c_int(Nvalues), DblArr(*times), DblArr(*values), c_int(interp))

    def interpolate(self, rebx, t):
        clibreboundx.rebx_interpolate.restype = c_double
        return clibreboundx.rebx_interpolate(byref(rebx), byref(self), c_double(t))

    def __del__(self):
        if self._b_needsfree_ == 1:
            clibreboundx.rebx_free_interpolator_pointers(byref(self))

Interpolator._fields_ = [  ("interpolation", c_int),
                    ("times", POINTER(c_double)),
                    ("values", POINTER(c_double)),
                    ("Nvalues", c_int),
                    ("y2", POINTER(c_double)),
                    ("klo", c_int)]

INTERPOLATION_TYPE = {"none":0, "spline":1}

# This list keeps pairing from C rebx_param_type enum to ctypes type 1-to-1. Derive the required mappings from it
REBX_C_TO_CTYPES = [["REBX_TYPE_NONE", None], ["REBX_TYPE_DOUBLE", c_double], ["REBX_TYPE_INT",c_int], ["REBX_TYPE_POINTER", c_void_p], ["REBX_TYPE_FORCE", Force], ["REBX_TYPE_UNIT32", c_uint32], ["REBX_TYPE_ORBIT", rebound.Orbit], ["REBX_TYPE_ODE", rebound.ODE], ["REBX_TYPE_VEC3D", rebound.Vec3d]]
REBX_CTYPES = {} # maps int value of rebx_param_type enum to ctypes type
REBX_C_PARAM_TYPES = {} # maps string of rebx_param_type enum to int
for i, pair in enumerate(REBX_C_TO_CTYPES):
    REBX_CTYPES[i] = pair[1]
    REBX_C_PARAM_TYPES[pair[0]] = i

from .params import Params
