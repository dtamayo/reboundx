from . import clibreboundx
from ctypes import Structure, c_double, POINTER, c_int, c_uint, c_long, c_ulong, c_void_p, c_char_p, CFUNCTYPE, byref, c_uint32, cast
import rebound
import reboundx


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

    #######################################
    # Functions for adding REBOUNDx effects
    #######################################

    def add(self, name):
        """
        This is the main function for adding effects to simulations.
        :param name: Name of the effect you wish to add.  See http://reboundx.readthedocs.io/en/latest/effects.html for a list of effects.
        :type name: string
        :rtype: rebx_effect structure
        """
        clibreboundx.rebx_add.restype = POINTER(rebx_effect)
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
        :rtype: rebx_effect structure
        """
        REBX_FUNCTION = CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(rebx_effect))
        self.custom_effects[function.__name__] = REBX_FUNCTION(function) # store a ref so it doesn't get garbage collected
        clibreboundx.rebx_add_custom_force.restype = POINTER(rebx_effect)
        ptr = clibreboundx.rebx_add_custom_force(byref(self), c_char_p(function.__name__.encode('ascii')), self.custom_effects[function.__name__], force_is_velocity_dependent)
        return ptr.contents

    def add_custom_post_timestep_modification(self, function):
        """
        This function allows you to add your own custom python function to REBOUNDx that is executed between integrator timesteps.
        You need to use this if you want to both use your own custom functions and the built-in REBOUNDx effects.  
        See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for details and a tutorial on how to use it.
        :param function: Custom Python function to be executed between timesteps.
        :type function: Function
        :type c: bool
        :rtype: rebx_effect structure 
        """
        REBX_FUNCTION = CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(rebx_effect))
        self.custom_effects[function.__name__] = REBX_FUNCTION(function) # store a ref so it doesn't get garbage collected
        clibreboundx.rebx_add_custom_post_timestep_modification.restype = POINTER(rebx_effect)
        ptr = clibreboundx.rebx_add_custom_post_timestep_modification(byref(self), c_char_p(function.__name__.encode('ascii')), self.custom_effects[function.__name__])
        return ptr.contents

    #######################################
    # Convenience Functions
    #######################################

    def rad_calc_beta(self, G, c, source_mass, source_luminosity, radius, density, Q_pr):
        """
        Calculates a particle's beta parameter (the ratio of the radiation force to the gravitational force).
        All values must be passed in the same units as used for the simulation as a whole (e.g., AU, Msun, yr/2pi).

        :param G: Gravitational constant
        :param c: Speed of light
        :param source_mass: Mass of radiation source
        :param source_luminosity: Luminosity of radiation source
        :param radius: grain's physical radius
        :param density: particle bulk density
        :param Q_pr: radiation pressure coefficient
        :type G: float
        :type c: float
        :type source_mass: float
        :type source_luminosity: float
        :type radius: float
        :type density: float
        :type Q_pr: float
        :rtype: float
        """
        clibreboundx.rebx_rad_calc_beta.restype = c_double
        return clibreboundx.rebx_rad_calc_beta(c_double(G), c_double(c), c_double(source_mass), c_double(source_luminosity), c_double(radius), c_double(density), c_double(Q_pr))

    def rad_calc_particle_radius(self, G, c, source_mass, source_luminosity, beta, density, Q_pr):
        """
        Calculates a particle's physical radius given its beta parameter.
        All values must be passed in the same units as used for the simulation as a whole (e.g., AU, Msun, yr/2pi).

        :param G: Gravitational constant
        :param c: Speed of light
        :param source_mass: Mass of radiation source
        :param source_luminosity: Luminosity of radiation source
        :param beta: Ratio of radiation force to gravitational force
        :param density: particle bulk density
        :param Q_pr: radiation pressure coefficient
        :type G: float
        :type c: float
        :type source_mass: float
        :type source_luminosity: float
        :type radius: float
        :type density: float
        :type Q_pr: float
        :rtype: float
        """
        clibreboundx.rebx_rad_calc_particle_radius.restype = c_double
        return clibreboundx.rebx_rad_calc_particle_radius(c_double(G), c_double(c), c_double(source_mass), c_double(source_luminosity), c_double(beta), c_double(density), c_double(Q_pr))
    
    def gr_hamiltonian(self, sim, params):
        """
        Calculates the value of the Hamiltonian from the particle states, when rebx_add_gr has been called.
        This also includes the classical Newtonian Hamiltonian.

        :param sim: REBOUND simulation.
        :param params: parameters instance returned by add_gr.
        :type sim: rebound.Simulation
        :type params: rebx_params_gr
        :rtype: float
        """
        clibreboundx.rebx_gr_hamiltonian.restype = c_double
        return clibreboundx.rebx_gr_hamiltonian(byref(sim), byref(params))
    
    def gr_potential_hamiltonian(self, sim, params):
        """
        Calculates the value of the Hamiltonian from the particle states, when rebx_add_gr has been called.
        This also includes the classical Newtonian Hamiltonian.

        :param sim: REBOUND simulation.
        :param params: parameters instance returned by add_gr_potential.
        :type sim: rebound.Simulation
        :type params: rebx_params_gr_potential
        :rtype: float
        """
        clibreboundx.rebx_gr_potential_hamiltonian.restype = c_double
        return clibreboundx.rebx_gr_potential_hamiltonian(byref(sim), byref(params))
    
    def gr_full_hamiltonian(self, sim, params):
        """
        Calculates the value of the Hamiltonian from the particle states, when rebx_add_gr_full has been called.
        This also includes the classical Newtonian Hamiltonian.

        :param sim: REBOUND simulation.
        :param params: parameters instance returned by add_gr_full.
        :type sim: rebound.Simulation
        :type params: rebx_params_gr_full
        :rtype: float
        """
        clibreboundx.rebx_gr_full_hamiltonian.restype = c_double
        return clibreboundx.rebx_gr_full_hamiltonian(byref(sim), byref(params))
    
    def tides_precession_hamiltonian(self, sim, params):
        clibreboundx.rebx_tides_precession_hamiltonian.restype = c_double
        return clibreboundx.rebx_tides_precession_hamiltonian(byref(sim), byref(params))

    
#################################################
# Generic REBOUNDx definitions
#################################################

class rebx_param(Structure): # need to define fields afterward because of circular ref in linked list
    pass    
rebx_param._fields_ =  [("paramPtr", c_void_p),
                        ("hash", c_uint32),
                        ("type_hash", c_uint32),
                        ("next", POINTER(rebx_param))]

class rebx_effect(Structure):
    @property 
    def params(self):
        params = Params(self)
        return params

rebx_effect._fields_ = [("object_type", c_uint32),
                        ("ap", POINTER(rebx_param)),
                        ("force", CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(rebx_effect))),
                        ("ptm", CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(rebx_effect))),
                        ("rebx", POINTER(Extras)),
                        ("next", POINTER(rebx_effect))]

class rebx_param_to_be_freed(Structure):
    pass
rebx_param_to_be_freed._fields_ =  [("param", POINTER(rebx_param)),
                                    ("next", POINTER(rebx_param_to_be_freed))]


# Need to put fields after class definition because of self-referencing
Extras._fields_ =  [("sim", POINTER(rebound.Simulation)),
                    ("effects", POINTER(rebx_effect)),
                    ("params_to_be_freed", POINTER(rebx_param_to_be_freed))]

from .params import Params
