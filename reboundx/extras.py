from . import clibreboundx
from ctypes import *
import rebound
c_default = 10064.915

class rebx_param(Structure): # need to define fields afterward because of circular ref in linked list
    pass    
rebx_param._fields_ = [("valPtr", c_void_p),
            ("param", c_int),
            ("next", POINTER(rebx_param))]

class rebx_param_to_be_freed(Structure):
    pass
rebx_param_to_be_freed._fields_ = [("param", POINTER(rebx_param)),
            ("next", POINTER(rebx_param_to_be_freed))]

class rebx_params_modify_orbits(Structure):
    _fields_ = [("p", c_double),
                ("coordinates", c_int)]

class rebx_params_gr(Structure):
    _fields_ = [("c", c_double)]

class rebx_params_radiation_forces(Structure):
    _fields_ = [("source", POINTER(rebound.Particle)),
                ("c", c_double)]

class Extras(Structure):
    """
    Main object used for all REBOUNDx operations, tied to a particular REBOUND simulation.
    This is an abstraction of the C struct rebx_extras, with all the C convenience functions
    and functions for adding effects implemented as methods of the class.  
    The fastest way to understand it is to follow the examples at :ref:`ipython_examples`.  
    """

    def __init__(self, sim):
        clibreboundx.rebx_initialize(byref(sim), byref(self)) # Use memory address ctypes allocated for rebx Structure in C
        self.add_Particle_props()
        self.coordinates = {"JACOBI":0, "BARYCENTRIC":1, "HELIOCENTRIC":2} # to use C version's REBX_COORDINATES enum
        sim._extras_ref = self # add a reference to this instance in sim to make sure it's not garbage collected

    def __del__(self):
        if self._b_needsfree_ == 1:
            clibreboundx.rebx_free_pointers(byref(self))

    def add_modify_orbits_direct(self):
        """
        Adds orbit modifications to the simulation, modifying orbital elements directly after each timestep.
        You must still set the relevant timescales on the individual particles.
        See :ref:`modules` for additional information, and for definitions of the relevant timescales.
        See :ref:`ipython_examples` for an example.
        """
        clibreboundx.rebx_add_modify_orbits_direct(byref(self))

    def add_modify_orbits_forces(self):
        """
        Adds orbit modifications to the simulation, implemented as forces that yield the desired effect.
        You must still set the relevant timescales on the individual particles.
        See :ref:`modules` for additional information, and for definitions of the relevant timescales.
        See :ref:`ipython_examples` for an example.
        """
        clibreboundx.rebx_add_modify_orbits_forces(byref(self))

    def check_c(self, c):
        if c is not None: # user passed c explicitly
            return c
      
        # c was not passed by user
         
        if self.sim.contents.G == 1: # if G = 1 (default) return default c
            return c_default
        else:
            raise ValueError("If you change G, you must pass c (speed of light) in appropriate units to add_gr, add_gr_potential, add_gr_full, and radiation_forces.  Setting the units in the simulation does not work with REBOUNDx.  See ipython_examples/GeneralRelativity.ipynb and ipython_examples/Radiation_Forces_Debris_Disk.ipynb")

    def add_gr(self, c=None):
        """
        Add general relativity corrections, treating only particles[0] as massive 
        (see :ref:`effectList` for details on the implementation). 
        Must pass the value of the speed of light if using non-default units (AU, Msun, yr/2pi)
        See :ref:`ipython_examples` for an example.
        """
        c = self.check_c(c)
        clibreboundx.rebx_add_gr(byref(self), c_double(c))
    
    def add_gr_full(self, c=None):
        """
        Add general relativity corrections, treating all particles as massive.
        (see :ref:`effectList` for details on the implementation). 
        Must pass the value of the speed of light if using non-default units (AU, Msun, yr/2pi)
        See :ref:`ipython_examples` for an example.
        """
        c = self.check_c(c)
        clibreboundx.rebx_add_gr_full(byref(self), c_double(c))

    def add_gr_potential(self, c=None):
        """
        Add general relativity corrections, using a simple potential that gets the precession right.
        (see :ref:`effectList` for details on the implementation). 
        Must pass the value of the speed of light if using non-default units (AU, Msun, yr/2pi)
        See :ref:`ipython_examples` for an example.
        """
        c = self.check_c(c)
        clibreboundx.rebx_add_gr_potential(byref(self), c_double(c))
    
    def add_radiation_forces(self, source, c=None):
        """
        Add radiation forces to the simulation (radiation pressure and Poynting-Robertson drag).
        (see :ref:`effectList` for details on the implementation). 
        Must pass the value of the speed of light if using non-default units (AU, Msun, yr/2pi),
        as well as the Particle in the Simulation that is the source of the radiation.
        See :ref:`ipython_examples` for an example.
        """
        c = self.check_c(c)
        clibreboundx.rebx_add_radiation_forces(byref(self), byref(source), c_double(c))

    def add_Particle_props(self):
        @property
        def tau_a(self):
            clibreboundx.rebx_get_tau_a.restype = c_double
            return clibreboundx.rebx_get_tau_a(byref(self))
        @tau_a.setter
        def tau_a(self, value):
            clibreboundx.rebx_set_tau_a(byref(self), c_double(value))
        @property
        def tau_e(self):
            clibreboundx.rebx_get_tau_e.restype = c_double
            return clibreboundx.rebx_get_tau_e(byref(self))
        @tau_e.setter
        def tau_e(self, value):
            clibreboundx.rebx_set_tau_e(byref(self), c_double(value))
        @property
        def tau_inc(self):
            clibreboundx.rebx_get_tau_inc.restype = c_double
            return clibreboundx.rebx_get_tau_inc(byref(self))
        @tau_inc.setter
        def tau_inc(self, value):
            clibreboundx.rebx_set_tau_inc(byref(self), c_double(value))
        @property
        def tau_omega(self):
            clibreboundx.rebx_get_tau_omega.restype = c_double
            return clibreboundx.rebx_get_tau_omega(byref(self))
        @tau_omega.setter
        def tau_omega(self, value):
            clibreboundx.rebx_set_tau_omega(byref(self), c_double(value))
        @property
        def tau_Omega(self):
            clibreboundx.rebx_get_tau_Omega.restype = c_double
            return clibreboundx.rebx_get_tau_Omega(byref(self))
        @tau_Omega.setter
        def tau_Omega(self, value):
            clibreboundx.rebx_set_tau_Omega(byref(self), c_double(value))
        @property
        def beta(self):
            clibreboundx.rebx_get_beta.restype = c_double
            return clibreboundx.rebx_get_beta(byref(self))
        @beta.setter
        def beta(self, value):
            clibreboundx.rebx_set_beta(byref(self), c_double(value))

        rebound.Particle.tau_a = tau_a
        rebound.Particle.tau_e = tau_e
        rebound.Particle.tau_inc = tau_inc
        rebound.Particle.tau_omega = tau_omega
        rebound.Particle.tau_Omega = tau_Omega
        rebound.Particle.beta = beta 

    def rad_calc_beta(self, particle_radius, density, Q_pr, L):
        """
        Calculates a particle's beta parameter (the ratio of the radiation force to the gravitational force) given
        the particle's physical radius, density, radiation pressure coefficient ``Q_pr``, and the star's luminosity.
        All values must be passed in the same units as used for the simulation as a whole (e.g., AU, Msun, yr/2pi).
        See the circumplanetary dust example in :ref:`ipython_examples`.
        """
        clibreboundx.rebx_rad_calc_beta.restype = c_double
        return clibreboundx.rebx_rad_calc_beta(byref(self), c_double(particle_radius), c_double(density), c_double(Q_pr), c_double(L))
    def rad_calc_particle_radius(self, beta, density, Q_pr, L):
        """
        Calculates a particle's physical radius given its beta parameter (the ratio of the radiation force to the gravitational force),
        density, radiation pressure coefficient ``Q_pr``, and the star's luminosity.
        All values must be passed in the same units as used for the simulation as a whole (e.g., AU, Msun, yr/2pi).
        """
        clibreboundx.rebx_rad_calc_particle_radius.restype = c_double
        return clibreboundx.rebx_rad_calc_particle_radius(byref(self), c_double(beta), c_double(density), c_double(Q_pr), c_double(L))


# Need to put fields after class definition because of self-referencing
Extras._fields_ = [("sim", POINTER(rebound.Simulation)),
                ("params_to_be_freed", POINTER(rebx_param_to_be_freed)),
                ("forces", POINTER(CFUNCTYPE(None, POINTER(rebound.Simulation)))),
                ("ptm", POINTER(CFUNCTYPE(None, POINTER(rebound.Simulation)))),
                ("Nptm", c_int),
                ("Nforces", c_int),
                ("modify_orbits_forces", rebx_params_modify_orbits),
                ("modify_orbits_direct", rebx_params_modify_orbits),
                ("gr", rebx_params_gr),
                ("radiation_forces", rebx_params_radiation_forces)]



