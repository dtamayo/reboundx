.. _modules:

REBx Effects & Parameters
=========================

.. _effectList:

List of REBOUNDx Effects
------------------------

Each of these can be added to a REBOUND simulation with the corresponding function in :ref:`add-effects`.  Each of these effects has at least one corresponding examples in both :ref:`c_examples` and :ref:`ipython_examples`. Make sure to look at the corresponding iPython example, which has further details on the physical model for the particular effect.

To see what particle parameters should be set for each effect, see :ref:`paramlist`.

*   *Modifying Orbital Elements*


    *   **modify_orbits_direct**
        

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

.. _gr:

gr
**

======================= ===============================================
Authors                 P. Shi, D. Tamayo, H. Rein
Implementation Paper    *In progress*
Based on                `Anderson et al. 1975 <http://labs.adsabs.harvard.edu/adsabs/abs/1975ApJ...200..221A/>`_.
C Example               :ref:`c_example_gr`
Python Example          `GeneralRelativity.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GeneralRelativity.ipynb>`_.
======================= ===============================================

This assumes that the masses are dominated by a single central body, and should be good enough for most applications with planets orbiting single stars.
It ignores terms that are smaller by of order the mass ratio with the central body.
It gets both the mean motion and precession correct, and will be significantly faster than :ref:`gr_full`, particularly with several bodies.

**Particle Parameters**

*None*

.. _gr_potential:

gr_potential
************

======================= ===============================================
Authors                 H. Rein, D. Tamayo
Implementation Paper    *In progress*
Based on                `Nobili and Roxburgh 1986 <http://labs.adsabs.harvard.edu/adsabs/abs/1986IAUS..114..105N/>`_.
C Example               :ref:`c_example_gr`
Python Example          `GeneralRelativity.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GeneralRelativity.ipynb>`_.
======================= ===============================================

This is the simplest potential you can use for general relativity.
It assumes that the masses are dominated by a single central body.
It gets the precession right, but gets the mean motion wrong by :math:`\mathcal{O}(GM/ac^2)`.  
It's the fastest option, and because it's not velocity-dependent, it automatically keeps WHFast symplectic.  
Nice if you have a single-star system, don't need to get GR exactly right, and want speed.

**Particle Parameters**

*None*

.. _gr_full:

gr_full
*******

======================= ===============================================
Authors                 P. Shi, H. Rein, D. Tamayo
Implementation Paper    *In progress*
Based on                `Newhall et al. 1983 <http://labs.adsabs.harvard.edu/adsabs/abs/1983A%26A...125..150N/>`_.
C Example               :ref:`c_example_gr`
Python Example          `GeneralRelativity.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GeneralRelativity.ipynb>`_.
======================= ===============================================

This algorithm incorporates GR effects from all bodies in the system, and is necessary for multiple massive bodies like stellar binaries.

**Particle Parameters**

*None*

Orbit Modifications
^^^^^^^^^^^^^^^^^^^

REBOUNDx offers two ways of modifying orbital elements (semimajor axis/eccentricity/inclination damping, precession, etc.)
In both cases, each particle is assigned evolution timescales for each orbital element.  
Positive timescales correspond to growth / progression, negative timescales correspond to damping / regression.  
Semimajor axes, eccentricities and inclinations grow / damp exponentially.  
Pericenters and nodes progress/regress linearly.

.. _modify_orbits_direct:

modify_orbits_direct
********************

======================= ===============================================
Authors                 D. Tamayo
Implementation Paper    *In progress*
Based on                `Lee & Peale 2002 <http://labs.adsabs.harvard.edu/adsabs/abs/2002ApJ...567..596L/>`_. 
C Example               :ref:`c_example_modify_orbits`
Python Example          `Migration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Migration.ipynb>`_
                        `EccAndIncDamping.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/EccAndIncDamping.ipynb>`_.
======================= ===============================================

This updates particles' positions and velocities between timesteps to achieve the desired changes to the osculating orbital elements (exponential growth/decay for a, e, inc, linear progression/regression for Omega/omega.
This nicely isolates changes to particular osculating elements, making it easier to interpret the resulting dynamics.  

**Particle Parameters**

======================= =========== ================================================
Parameter name          Required    Description
======================= =========== ================================================
tau_a                   No          Semimajor axis exponential growth/damping timescale
tau_e                   No          Eccentricity exponential growth/damping timescale
tau_inc                 No          Inclination axis exponential growth/damping timescale
tau_Omega               No          Period of linear nodal precession/regression
tau_omega               No          Period of linear apsidal precession/regression
======================= =========== ================================================

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

