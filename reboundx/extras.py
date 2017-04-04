from . import clibreboundx
from ctypes import Structure, c_double, POINTER, c_int, c_uint, c_long, c_ulong, c_void_p, c_char_p, CFUNCTYPE, byref, c_uint32, c_uint, cast, c_char
import rebound
import reboundx
import warnings

INTEGRATORS = {"implicit_midpoint": 0, "rk4":1, "euler": 2, "none": 3}
  
REBX_BINARY_WARNINGS = [
        ("REBOUNDx Error: Cannot read binary file. Check filename and file contents.", 1),
        ("REBOUNDx Error: Binary file was corrupt. Could not read.", 2),
        ("REBOUNDx Warning: Binary file was saved with a different version of REBOUNDx. Binary format might have changed.", 4),
        ("REBOUNDx Warning: At least one parameter in the binary file was not loaded. Check simulation.", 8),
        ("REBOUNDx Warning: At least one particle's parameters in the binary file were not loaded. Check simulation.", 16),
        ("REBOUNDx Warning: At least one effect and its parameters were not loaded from the binary file. Check simulation.", 32),
        ("REBOUNDx Warning: At least one field in the binary field was not recognized, and not loaded. Probably binary was created with more recent REBOUNDx version than you are using.", 64),
]

class Extras(Structure):
    """
    Main object used for all REBOUNDx operations, tied to a particular REBOUND simulation.
    This is an abstraction of the C struct rebx_extras, with all the C convenience functions
    and functions for adding effects implemented as methods of the class.  
    The fastest way to understand it is to follow the examples at :ref:`ipython_examples`.  
    """

    def __init__(self, sim):
        #first check whether additional_forces or post_timestep_modifications is set on sim.  If so, raise error
        #if cast(sim._additional_forces, c_void_p).value is not None or cast(sim._post_timestep_modifications, c_void_p).value is not None:
        #    raise AttributeError("sim.additional_forces or sim.post_timestep_modifications was already set.  If you want to use REBOUNDx, you need to add custom effects through REBOUNDx.  See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for a tutorial.")
        
        clibreboundx.rebx_initialize(byref(sim), byref(self)) # Use memory address ctypes allocated for rebx Structure in C
        if not hasattr(sim, "_extras_ref"): # if REBOUNDx wasn't already attached, check for warnings in case additional_forces or ptm were already set.
            sim.process_messages()
        sim._extras_ref = self # add a reference to this instance in sim to make sure it's not garbage collected
        self.custom_effects = {} # dictionary to keep references to custom effects so they don't get garbage collected

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
            clibreboundx.rebx_free_effects(byref(self))
            clibreboundx.rebx_free_params(byref(self))

    def remove_from_simulation(self, sim):
        """
        Disattaches reboundx from simulation, removing all effects.
        """
        del sim._extras_ref
        clibreboundx.rebx_remove_from_simulation(byref(sim))

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

    def add(self, name):
        """
        This is the main function for adding effects to simulations.
        :param name: Name of the effect you wish to add.  See http://reboundx.readthedocs.io/en/latest/effects.html for a list of effects.
        :type name: string
        :rtype: Effect structure
        """
        clibreboundx.rebx_add.restype = POINTER(Effect)
        ptr = clibreboundx.rebx_add(byref(self), c_char_p(name.encode('ascii')))
        self.sim.contents.process_messages()
        return ptr.contents

    def add_custom_force(self, function, force_is_velocity_dependent):
        """
        This function allows you to add your own custom python function to REBOUNDx that updates particle accelerations.  
        You need to use this if you want to both use your own custom functions and the built-in REBOUNDx effects.  
        See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for details and a tutorial on how to use it.
        :param function: Custom Python function for updating particle accelerations.
        :param force_is_velocity_dependent: Whether your custom force depends on the particle velocities.
        :type function: Function
        :type c: bool
        :rtype: Effect structure
        """
        REBX_FORCE = CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(Effect), POINTER(rebound.Particle), c_int)
        self.custom_effects[function.__name__] = REBX_FORCE(function) # store a ref so it doesn't get garbage collected
        clibreboundx.rebx_add_custom_force.restype = POINTER(Effect)
        ptr = clibreboundx.rebx_add_custom_force(byref(self), c_char_p(function.__name__.encode('ascii')), self.custom_effects[function.__name__], force_is_velocity_dependent)
        return ptr.contents

    def add_custom_operator(self, function):
        """
        This function allows you to add your own custom python function to REBOUNDx that is executed between integrator timesteps.
        You need to use this if you want to both use your own custom functions and the built-in REBOUNDx effects.  
        See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for details and a tutorial on how to use it.
        :param function: Custom Python function to be executed between timesteps.
        :type function: Function
        :type c: bool
        :rtype: Effect structure 
        """
        REBX_OPERATOR = CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(Effect), c_double, c_int)
        self.custom_effects[function.__name__] = REBX_OPERATOR(function) # store a ref so it doesn't get garbage collected
        clibreboundx.rebx_add_custom_operator.restype = POINTER(Effect)
        ptr = clibreboundx.rebx_add_custom_operator(byref(self), c_char_p(function.__name__.encode('ascii')), self.custom_effects[function.__name__])
        return ptr.contents
    
    def get_effect(self, name):
        clibreboundx.rebx_get_effect.restype = POINTER(Effect)
        ptr = clibreboundx.rebx_get_effect(byref(self), c_char_p(name.encode('ascii')))
        if ptr:
            return ptr.contents
        else:
            warnings.warn("Parameter {0} not found".format(name), RuntimeWarning)
            return

    def gr_acc(self, acc, C2):
        #sixdoub = c_double*6
        #clibreboundx.rebx_gr_acc.restype = POINTER(sixdoub)
        clibreboundx.rebx_gr_acc(byref(self), acc, c_double(C2))
    def calculate_energy(self):
        clibreboundx.rebx_calculate_energy.restype = c_double
        return clibreboundx.rebx_calculate_energy(self.sim)
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
        self.sim.contents.process_messages()
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
    
#################################################
# Generic REBOUNDx definitions
#################################################

class Param(Structure): # need to define fields afterward because of circular ref in linked list
    pass    
Param._fields_ =  [("contents", c_void_p),
                        ("hash", c_uint32),
                        ("name", c_char_p),
                        ("param_type", c_int),
                        ("python_type", c_int),
                        ("ndim", c_int),
                        ("shape", POINTER(c_int)),
                        ("strides", POINTER(c_int)),
                        ("size", c_int),
                        ("next", POINTER(Param))]

class Effect(Structure):
    @property 
    def params(self):
        params = Params(self)
        return params

Effect._fields_ = [ ("hash", c_uint32),
                    ("name", c_char_p),
                    ("ap", POINTER(Param)),
                    ("force_as_operator", c_int),
                    ("force", CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(Effect), POINTER(rebound.Particle), c_int)),
                    ("operator", CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(Effect), c_double, c_int)),
                    ("operator_order", c_int),
                    ("rebx", POINTER(Extras)),
                    ("next", POINTER(Effect)),
                    ("prev", POINTER(Effect)),
                    ("pad", c_char*100)]

class Param_to_be_freed(Structure):
    pass
Param_to_be_freed._fields_ =  [("param", POINTER(Param)),
                                    ("next", POINTER(Param_to_be_freed))]


# Need to put fields after class definition because of self-referencing
Extras._fields_ =  [("sim", POINTER(rebound.Simulation)),
                    ("effects", POINTER(Effect)),
                    ("params_to_be_freed", POINTER(Param_to_be_freed)),
                    ("_integrator", c_int),]

from .params import Params
