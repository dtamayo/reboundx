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

        for i in range(sim.N):
            sim.particles[i].ap = None

    def __del__(self):
        if self._b_needsfree_ == 1:
            clibreboundx.rebx_free_effects(byref(self))
            clibreboundx.rebx_free_params(byref(self))
    
    #######################################
    # Functions for adding REBOUNDx effects
    #######################################

    def add_effect(self, name):
        clibreboundx.rebx_add_effect.restype = POINTER(rebx_effect)
        ptr = clibreboundx.rebx_add_effect(byref(self), c_char_p(name.encode('ascii')))
        return ptr.contents
    
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

