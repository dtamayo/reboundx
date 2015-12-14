.. _modules:

REBx Effects & Parameters
=========================

.. _effectList:

List of REBOUNDx Effects
------------------------

Each of these can be added to a REBOUND simulation with the corresponding function in :ref:`add-effects`.  Each of these effects has at least one corresponding examples in both :ref:`c_examples` and :ref:`ipython_examples`. Make sure to look at the corresponding iPython example, which has further details on the physical model for the particular effect.

To see what particle parameters should be set for each effect, see :ref:`paramlist`.

* **modify_orbits_direct**
  
  This updates particles' positions and velocities between timesteps to achieve the desired changes to the osculating orbital elements.  This is the approach of Lee & Peale (2002).  This nicely isolates changes to particular osculating elements, making it easier to interpret the resulting dynamics.  

  Each particle is assigned evolution timescales for each orbital element.  Positive timescales correspond to growth / progression, negative timescales correspond to damping / regression.  Semimajor axes, eccentricities and inclinations grow / damp exponentially.  Pericenters and nodes progress/regress linearly.

* **modify_orbits_forces**

  This adds additional forces that, averaged over an orbit, yield exponential growth/damping of the semimajor axes, eccentricities and inclinations (set positive timescale for growth, negative for damping).  These equations follow Papaloizou & Larwood (2000).  These provide more physically realistic coupling of the orbital elements' evolutions to one another, e.g., adding an eccentricity damping timescale causes semimajor axis evolution at order :math:`e^2` (equivalent to a p parameter of 1 for modify_orbit_direct), cf. Goldreich & Schlichting (2014), Deck & Batygin (2015), Delisle et al. (2015). 

* **gr**

  Adds post-Newtonian corrections to a simulation, treating only ``particles[0]`` as massive for the modification.  This should be OK for objects orbiting a single star.  Follows Eq. 2 of Benitez and Gallardo 2008.  It gets both the mean motion and precession correct, and will be significantly faster than gr_full, particularly with several bodies.

* **gr_full**

  Follows Eq. 1 of Benitez and Gallardo 2008, treating all objects in the simulation as massive for the post-Newtonian correction.  This is the most general / safest implementation, but will be slower.

* **gr_potential**

  This is the simplest potential you can use for GR (Nobili & Roxburgh 1986). It gets the precession right, but gets the mean motion wrong by :math:`\mathcal{O}(GM/ac^2)`.  It's the fasest option, and because it's not velocity-dependent, it automatically keeps WHFast symplectic.  Nice if you don't need to get GR exactly right and want speed.

* **radiation_forces**

  This adds radiation forces, following Eq. 5 of Burns, Lamy & Soter (1979).  This incorporates both radiation pressure and Poynting-Robertson drag.  

.. _paramlist:

List of Particle Parameters
---------------------------

=============== ========================================= ============================================ 
Parameter name  Effect                                    Description
=============== ========================================= ============================================ 
tau_a           modify_orbit_direct, modify_orbits_forces Semimajor axis exponential growth/damping timescale
tau_e           modify_orbit_direct, modify_orbits_forces Eccentricity exponential growth/damping timescale
tau_inc         modify_orbit_direct, modify_orbits_forces Inclination axis exponential growth/damping timescale
tau_Omega       modify_orbit_direct, modify_orbits_forces Period of linear nodal precession/regression
tau_omega       modify_orbit_direct, modify_orbits_forces Period of linear apsidal precession/regression
beta            radiation_forces                          Ratio of radiation to gravitational force (Burns et al. 1979)
=============== ========================================= ============================================ 

