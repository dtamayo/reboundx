from . import clibreboundx
from ctypes import Structure, c_double, POINTER, c_int, c_uint, c_long, c_ulong, c_void_p, c_char_p, CFUNCTYPE, byref, c_uint32, cast
import rebound
import reboundx

C_DEFAULT = 10064.915 # speed of light in default units (G = 1) of AU / (yr/2pi)

class Extras(Structure):
    """
    Main object used for all REBOUNDx operations, tied to a particular REBOUND simulation.
    This is an abstraction of the C struct rebx_extras, with all the C convenience functions
    and functions for adding effects implemented as methods of the class.  
    The fastest way to understand it is to follow the examples at :ref:`ipython_examples`.  
    """

    def __init__(self, sim):
        #first check whether additional_forces or post_timestep_modifications is set on sim.  If so, raise error
        if cast(sim._additional_forces, c_void_p).value is not None or cast(sim._post_timestep_modifications, c_void_p).value is not None:
            raise AttributeError("sim.additional_forces or sim.post_timestep_modifications was already set.  If you want to use REBOUNDx, you need to add custom effects through REBOUNDx.  See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for a tutorial.")
        
        clibreboundx.rebx_initialize(byref(sim), byref(self)) # Use memory address ctypes allocated for rebx Structure in C
        sim._extras_ref = self # add a reference to this instance in sim to make sure it's not garbage collected
        self.custom_effects = {} # dictionary to keep references to custom effects so they don't get garbage collected

    def __del__(self):
        if self._b_needsfree_ == 1:
            clibreboundx.rebx_free_effects(byref(self))
            clibreboundx.rebx_free_params(byref(self))
    
    #######################################
    # Functions for adding REBOUNDx effects
    #######################################

    def add_modify_orbits_direct(self):
        """
        :rtype: rebx_params_modify_orbits_direct
        """
        clibreboundx.rebx_add_modify_orbits_direct.restype = POINTER(rebx_params_modify_orbits_direct)
        return clibreboundx.rebx_add_modify_orbits_direct(byref(self)).contents

    def add_modify_orbits_forces(self):
        """
        :rtype: rebx_params_modify_orbits_forces
        """
        clibreboundx.rebx_add_modify_orbits_forces.restype = POINTER(rebx_params_modify_orbits_forces)
        return clibreboundx.rebx_add_modify_orbits_forces(byref(self)).contents

    def add_gr(self, source_index=0, c=C_DEFAULT): 
        """
        You must pass c (the speed of light) in whatever units you choose if you don't use default units of AU, (yr/2pi) and Msun.
        
        :param source_index: Index in the particles array of the massive body that is the source of the GR corrections.
        :param c: Speed of light in appropriate units.
        :type source_index: int
        :type c: float
        :rtype: rebx_params_gr
        """
        clibreboundx.rebx_add_gr.restype = POINTER(rebx_params_gr)
        return clibreboundx.rebx_add_gr(byref(self), c_int(source_index), c_double(c)).contents
    
    def add_gr_full(self, c=C_DEFAULT):
        """
        You must pass c (the speed of light) in whatever units you choose if you don't use default units of AU, (yr/2pi) and Msun.
        
        :param c: Speed of light in appropriate units.
        :type c: float
        :rtype: rebx_params_gr_full
        """
        clibreboundx.rebx_add_gr_full.restype = POINTER(rebx_params_gr_full)
        return clibreboundx.rebx_add_gr_full(byref(self), c_double(c)).contents

    def add_gr_potential(self, source_index=0, c=C_DEFAULT):
        """
        You must pass c (the speed of light) in whatever units you choose if you don't use default units of AU, (yr/2pi) and Msun.
        
        
        :param source_index: Index in the particles array of the massive body that is the source of the GR corrections.
        :param c: Speed of light in appropriate units.
        :type source_index: int
        :type c: float
        :rtype: rebx_params_gr_potential
        """
        clibreboundx.rebx_add_gr_potential.restype = POINTER(rebx_params_gr_potential)
        return clibreboundx.rebx_add_gr_potential(byref(self), c_int(source_index), c_double(c)).contents
    
    def add_radiation_forces(self, source_index=0, c=C_DEFAULT):
        """
        You must pass c (the speed of light) in whatever units you choose if you don't use default units of AU, (yr/2pi) and Msun.
        
        :param source_index: Index in the particles array of the body that is the source of the radiation.
        :param c: Speed of light in appropriate units.
        :type source_index: int
        :type c: float
        :rtype: rebx_params_radiation_forces
        """
        clibreboundx.rebx_add_radiation_forces.restype = POINTER(rebx_params_radiation_forces)
        return clibreboundx.rebx_add_radiation_forces(byref(self), c_int(source_index), c_double(c)).contents
    
    def add_modify_mass(self):
        """
        :rtype: None
        """
        clibreboundx.rebx_add_modify_mass(byref(self))

    def add_custom_force(self, function, force_is_velocity_dependent):
        """
        This function allows you to add your own custom python function to REBOUNDx that updates particle accelerations.  
        You need to use this if you want to both use your own custom functions and the built-in REBOUNDx effects.  
        See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for details and a tutorial on how to use it.

        :param function: Custom Python function for updating particle accelerations.
        :param force_is_velocity_dependent: Whether your custom force depends on the particle velocities.
        :type function: Function
        :type c: bool
        :rtype: None
        """
        REBX_FUNCTION = CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(rebx_effect))
        self.custom_effects[function.__name__] = REBX_FUNCTION(function) # store a ref so it doesn't get garbage collected
        clibreboundx.rebx_add_custom_force(byref(self), self.custom_effects[function.__name__], force_is_velocity_dependent, None) # set custom_params ptr to NULL (not supported in python version)
    
    def add_custom_post_timestep_modification(self, function):
        """
        This function allows you to add your own custom python function to REBOUNDx that is executed between integrator timesteps.
        You need to use this if you want to both use your own custom functions and the built-in REBOUNDx effects.  
        See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for details and a tutorial on how to use it.

        :param function: Custom Python function to be executed between timesteps.
        :type function: Function
        :rtype: None
        """
        REBX_FUNCTION = CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(rebx_effect))
        self.custom_effects[function.__name__] = REBX_FUNCTION(function) # store a ref so it doesn't get garbage collected
        clibreboundx.rebx_add_custom_post_timestep_modification(byref(self), self.custom_effects[function.__name__], None) # set custom_params ptr to NULL (not supported in python version)

    #######################################
    # Convenience Functions
    #######################################

    def rad_calc_beta(self, params, particle_radius, density, Q_pr, L):
        """
        Calculates a particle's beta parameter (the ratio of the radiation force to the gravitational force).
        All values must be passed in the same units as used for the simulation as a whole (e.g., AU, Msun, yr/2pi).

        :param params: parameters instance returned by add_radiation_forces.
        :param particle_radius: grain's physical radius
        :param density: particle bulk density
        :param Q_pr: radiation pressure coefficient
        :param L: Radiation source's luminosity
        :type params: rebx_params_radiation_forces
        :type particle_radius: float
        :type density: float
        :type Q_pr: float
        :type L: float
        :rtype: float
        """
        clibreboundx.rebx_rad_calc_beta.restype = c_double
        return clibreboundx.rebx_rad_calc_beta(byref(self), byref(params), c_double(particle_radius), c_double(density), c_double(Q_pr), c_double(L))

    def rad_calc_particle_radius(self, params, beta, density, Q_pr, L):
        """
        Calculates a particle's physical radius given its beta parameter.
        All values must be passed in the same units as used for the simulation as a whole (e.g., AU, Msun, yr/2pi).

        :param params: parameters instance returned by add_radiation_forces.
        :param beta: ratio of radiation pressure force to gravitational force from the radiation source.
        :param density: particle bulk density
        :param Q_pr: radiation pressure coefficient
        :param L: Radiation source's luminosity
        :type params: rebx_params_radiation_forces
        :type beta: float
        :type density: float
        :type Q_pr: float
        :type L: float
        :rtype: float
        """
        clibreboundx.rebx_rad_calc_particle_radius.restype = c_double
        return clibreboundx.rebx_rad_calc_particle_radius(byref(self), byref(params), c_double(beta), c_double(density), c_double(Q_pr), c_double(L))
    
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
    
#######################################
# Effect parameter class definitions
#######################################

class rebx_params_modify_orbits_direct(Structure):
    _fields_ = [("p", c_double),
                ("coordinates", c_int)]

class rebx_params_modify_orbits_forces(Structure):
    _fields_ = [("coordinates", c_int)]

class rebx_params_gr(Structure):
    _fields_ = [("source_index", c_int),
                ("c", c_double)]

class rebx_params_gr_potential(Structure):
    _fields_ = [("source_index", c_int),
                ("c", c_double)]

class rebx_params_gr_full(Structure):
    _fields_ = [("c", c_double)]

class rebx_params_radiation_forces(Structure):
    _fields_ = [("source_index", c_int),
                ("c", c_double)]

#################################################
# Generic REBOUNDx definitions
#################################################

class rebx_param(Structure): # need to define fields afterward because of circular ref in linked list
    pass    
rebx_param._fields_ = [("paramPtr", c_void_p),
            ("param_type", c_int),
            ("next", POINTER(rebx_param))]

class rebx_effect(Structure):
    pass
rebx_effect._fields_ = [("paramsPtr", c_void_p),
                        ("functionPtr", CFUNCTYPE(None, POINTER(rebound.Simulation), POINTER(rebx_effect))),
                        ("effect_type", c_uint32),
                        ("is_force", c_int),
                        ("next", POINTER(rebx_effect))]

class rebx_param_to_be_freed(Structure):
    pass
rebx_param_to_be_freed._fields_ = [("param", POINTER(rebx_param)),
            ("next", POINTER(rebx_param_to_be_freed))]


# Need to put fields after class definition because of self-referencing
Extras._fields_ = [("sim", POINTER(rebound.Simulation)),
                ("effects", POINTER(rebx_effect)),
                ("params_to_be_freed", POINTER(rebx_param_to_be_freed))]
