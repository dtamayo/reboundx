.. _modules:

REBx Effects & Parameters
=========================

.. _effectList:

List of REBOUNDx Effects
------------------------

Each of these can be added to a REBOUND simulation with the corresponding function in :ref:`add-effects`.  Each of these effects has at least one corresponding examples in both :ref:`c_examples` and :ref:`ipython_examples`. Make sure to look at the corresponding iPython example, which has further details on the physical model for the particular effect.

To see what particle parameters should be set for each effect, see :ref:`paramlist`.

*   *Modifying Orbital Elements*

    REBOUNDx offers two ways of modifying orbital elements (semimajor axis/eccentricity/inclination damping, precession, etc.)
    In both cases, each particle is assigned evolution timescales for each orbital element.  
    Positive timescales correspond to growth / progression, negative timescales correspond to damping / regression.  
    Semimajor axes, eccentricities and inclinations grow / damp exponentially.  
    Pericenters and nodes progress/regress linearly.

    *   **modify_orbits_direct**
        
        This updates particles' positions and velocities between timesteps to achieve the desired changes to the osculating orbital elements.  
        This is the approach of Lee & Peale (2002).  
        This nicely isolates changes to particular osculating elements, making it easier to interpret the resulting dynamics.  

    *   **modify_orbits_forces**
        
        This adds additional forces that, averaged over an orbit, yield exponential growth/damping of the 
        semimajor axes, eccentricities and inclinations (set positive timescale for growth, negative for damping).  
        These equations follow Papaloizou & Larwood (2000).  
        These provide more physically realistic coupling of the orbital elements' evolutions to one another, 
        e.g., adding an eccentricity damping timescale causes semimajor axis evolution at order :math:`e^2` 
        (equivalent to a p parameter of 1 for modify_orbit_direct), 
        cf. Goldreich & Schlichting (2014), Deck & Batygin (2015), Delisle et al. (2015). 
        Eccentricity damping also induces pericenter precession.

*   *General Relativity Corrections*

    *   **gr**

        Adds post-Newtonian corrections to a simulation, treating only ``particles[0]`` as massive for the modification.  
        This should be OK for objects orbiting a single star.  
        Follows Eq. 2 of Benitez and Gallardo 2008.  
        It gets both the mean motion and precession correct, and will be significantly faster than gr_full, particularly with several bodies.

    *   **gr_full**

        Follows Eq. 1 of Benitez and Gallardo 2008, treating all objects in the simulation as massive for the post-Newtonian correction.  
        This is the most general / safest implementation, but will be slower.

    *   **gr_potential**

        This is the simplest potential you can use for GR (Nobili & Roxburgh 1986). 
        It gets the precession right, but gets the mean motion wrong by :math:`\mathcal{O}(GM/ac^2)`.  
        It's the fasest option, and because it's not velocity-dependent, it automatically keeps WHFast symplectic.  
        Nice if you don't need to get GR exactly right and want speed.

*   **radiation_forces**

    This adds radiation forces, following Eq. 5 of Burns, Lamy & Soter (1979).  
    It incorporates both radiation pressure and Poynting-Robertson drag.  

.. _effectDescriptions:

Effect Descriptions
-------------------

General Relativity
^^^^^^^^^^^^^^^^^^

gr
**

gr_potential
************

+-----------------------+-----------------------------------------------+
| Authors               | H. Rein, D. Tamayo                            |
+=======================+===============================================+
| Authors               | H. Rein, D. Tamayo                            |
+-----------------------+-----------------------------------------------+
| Authors               | H. Rein, D. Tamayo                            |
+-----------------------+-----------------------------------------------+

===================== ===============================================
Authors               H. Rein, D. Tamayo
===================== ===============================================
Implementation Paper  *In progress*
===================== ===============================================

**Authors**: H. Rein, D. Tamayo

**Implementation paper**: *In progress*

**Based on**: `Nobili and Roxburgh 1986 <http://labs.adsabs.harvard.edu/adsabs/abs/1986IAUS..114..105N/>`_.

**C Example**:
        
**Python Example**: `GeneralRelativity.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GeneralRelativity.ipynb>`_.

This is the simplest potential you can use for general relativity.
It adds corrections from a single massive body, specified by `source`.
It gets the precession right, but gets the mean motion wrong by :math:`\mathcal{O}(GM/ac^2)`.  
It's the fastest option, and because it's not velocity-dependent, it automatically keeps WHFast symplectic.  
Nice if you don't need to get GR exactly right and want speed.

**Particle Parameters**

*None*

===================== ======== ===============================================
Name                  Required Description
===================== ======== ===============================================
===================== ======== ===============================================

gr_full
*******


.. _paramlist:

Table of Particle Parameters
----------------------------

In the Python version, particle parameters can be accessed with ``sim.particles[1].param`` or set with ``sim.particles[1].param = value``, where ``param`` is a parameter name from the table below.  

In the C version, you get with ``rebx_get_param(&sim->particles[3])`` and set with ``rebx_set_param(&sim->particles[3], value)``.

There are plenty of examples in :ref:`c_examples` and :ref:`ipython_examples`.

=============== ========================================= ============================================ 
Parameter name  Effect                                    Description
=============== ========================================= ============================================ 
tau_a           modify_orbit_direct, modify_orbits_forces Semimajor axis exponential growth/damping timescale
tau_e           modify_orbit_direct, modify_orbits_forces Eccentricity exponential growth/damping timescale
tau_inc         modify_orbit_direct, modify_orbits_forces Inclination axis exponential growth/damping timescale
tau_Omega       modify_orbit_direct                       Period of linear nodal precession/regression
tau_omega       modify_orbit_direct                       Period of linear apsidal precession/regression
beta            radiation_forces                          Ratio of radiation to gravitational force (Burns et al. 1979)
=============== ========================================= ============================================ 

