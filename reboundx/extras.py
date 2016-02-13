from . import clibreboundx
from ctypes import Structure, c_double, POINTER, c_int, c_uint, c_long, c_ulong, c_void_p, c_char_p, CFUNCTYPE, byref, c_uint32
from .tools import check_c
import rebound

class Extras(Structure):
    """
    Main object used for all REBOUNDx operations, tied to a particular REBOUND simulation.
    This is an abstraction of the C struct rebx_extras, with all the C convenience functions
    and functions for adding effects implemented as methods of the class.  
    The fastest way to understand it is to follow the examples at :ref:`ipython_examples`.  
    """

    def __init__(self, sim):
        clibreboundx.rebx_initialize(byref(sim), byref(self)) # Use memory address ctypes allocated for rebx Structure in C
        sim._extras_ref = self # add a reference to this instance in sim to make sure it's not garbage collected

    def __del__(self):
        if self._b_needsfree_ == 1:
            clibreboundx.rebx_free_pointers(byref(self))
    
    def add_Particle_props(self):
        print("monkey")

    #######################################
    # Functions for adding REBOUNDx effects
    #######################################

    def add_modify_orbits_direct(self):
        """
        Adds orbit modifications to the simulation, modifying orbital elements directly after each timestep.
        You must still set the relevant timescales on the individual particles.
        See :ref:`modules` for additional information, and for definitions of the relevant timescales.
        See :ref:`ipython_examples` for an example.
        """
        clibreboundx.rebx_add_modify_orbits_direct.restype = POINTER(rebx_params_modify_orbits_direct)
        return clibreboundx.rebx_add_modify_orbits_direct(byref(self)).contents

    def add_modify_orbits_forces(self):
        """
        Adds orbit modifications to the simulation, implemented as forces that yield the desired effect.
        You must still set the relevant timescales on the individual particles.
        See :ref:`modules` for additional information, and for definitions of the relevant timescales.
        See :ref:`ipython_examples` for an example.
        """
        clibreboundx.rebx_add_modify_orbits_forces.restype = POINTER(rebx_params_modify_orbits_forces)
        return clibreboundx.rebx_add_modify_orbits_forces(byref(self)).contents

    def add_gr(self, source=None, c=None):
        """
        Add general relativity corrections from a single body, specified by source (defaults to particles[0]).
        (see :ref:`effectList` for details on the implementation). 
        Must pass the value of the speed of light if using non-default units (AU, Msun, yr/2pi)
        See :ref:`ipython_examples` for an example.
        """
        c = check_c(self, c)
        if source is not None:
            source = byref(source)
        clibreboundx.rebx_add_gr.restype = POINTER(rebx_params_gr)
        return clibreboundx.rebx_add_gr(byref(self), source, c_double(c)).contents # Sets source to particles[0] in C code when passed NULL (=None)
    
    def add_gr_full(self, c=None):
        """
        Add general relativity corrections, treating all particles as massive.
        (see :ref:`effectList` for details on the implementation). 
        Must pass the value of the speed of light if using non-default units (AU, Msun, yr/2pi)
        See :ref:`ipython_examples` for an example.
        """
        c = check_c(self, c)
        clibreboundx.rebx_add_gr_full.restype = POINTER(rebx_params_gr_full)
        return clibreboundx.rebx_add_gr_full(byref(self), c_double(c)).contents

    def add_gr_potential(self, source=None, c=None):
        """
        
        :param source: Source of GR effect. Defaults to sim.particles[0].
        :param c: Speed of light in appropriate units. Defaults to AU/(yr/2pi)
        :type source: rebound.Particle
        :type c: double
        :rtype: rebx_params_gr_potential

        """
        c = check_c(self, c)
        if source is not None:
            source = byref(source)
        clibreboundx.rebx_add_gr_potential.restype = POINTER(rebx_params_gr_potential)
        return clibreboundx.rebx_add_gr_potential(byref(self), source, c_double(c)).contents
    
    def add_radiation_forces(self, source=None, c=None):
        """
        Add radiation forces to the simulation (radiation pressure and Poynting-Robertson drag).
        (see :ref:`effectList` for details on the implementation). 
        Must pass the value of the speed of light if using non-default units (AU, Msun, yr/2pi),
        as well as the Particle in the Simulation that is the source of the radiation.
        See :ref:`ipython_examples` for an example.
        """
        c = check_c(self, c)
        if source is not None:
            source = byref(source)
        clibreboundx.rebx_add_radiation_forces.restype = POINTER(rebx_params_radiation_forces)
        return clibreboundx.rebx_add_radiation_forces(byref(self), source, c_double(c)).contents
    
    #######################################
    # Convenience Functions
    #######################################

    def rad_calc_beta(self, params, particle_radius, density, Q_pr, L):

        """
        Calculates a particle's beta parameter (the ratio of the radiation force to the gravitational force) given
        the params returned from add_radiation_forces, the particle's physical radius, density, radiation pressure 
        coefficient ``Q_pr``, and the star's luminosity.
        All values must be passed in the same units as used for the simulation as a whole (e.g., AU, Msun, yr/2pi).
        See the circumplanetary dust example in :ref:`ipython_examples`.
        """
        clibreboundx.rebx_rad_calc_beta.restype = c_double
        return clibreboundx.rebx_rad_calc_beta(byref(self), byref(params), c_double(particle_radius), c_double(density), c_double(Q_pr), c_double(L))
    def rad_calc_particle_radius(self, params, beta, density, Q_pr, L):
        """
        Calculates a particle's physical radius given the params returned from add_radiation_forces, the particle's
        beta parameter (the ratio of the radiation force to the gravitational force),
        density, radiation pressure coefficient ``Q_pr``, and the star's luminosity.
        All values must be passed in the same units as used for the simulation as a whole (e.g., AU, Msun, yr/2pi).
        """
        clibreboundx.rebx_rad_calc_particle_radius.restype = c_double
        return clibreboundx.rebx_rad_calc_particle_radius(byref(self), byref(params), c_double(beta), c_double(density), c_double(Q_pr), c_double(L))
    
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
    _fields_ = [("source_index", c_int), # testing
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
