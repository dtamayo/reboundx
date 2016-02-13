.. _effects:

REBx Effects & Parameters
=========================

Each of these effects can be added to a REBOUND simulation with the corresponding function in :ref:`add-effects`.  Each of these effects has at least one corresponding examples in both :ref:`c_examples` and :ref:`ipython_examples`. Make sure to look at the corresponding iPython example, which has further details on the physical model for the particular effect.

In the Python version, particle parameters can be accessed with ``sim.particles[1].param`` or set with ``sim.particles[1].param = value``, where ``param`` is a parameter name from the table below.  

In the C version, you get with 

There are plenty of examples in :ref:`c_examples` and :ref:`ipython_examples`.
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

**Effect Structure**: *rebx_params_gr*

=========================== ==================================================================
Field (C type)              Description
=========================== ==================================================================
c (double)                  Speed of light in the units used for the simulation.
source_index (int)          Index in the `particles` array for the massive central body.
=========================== ==================================================================

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

**Effect Structure**: *rebx_params_gr_potential*

=========================== ==================================================================
Field (C type)              Description
=========================== ==================================================================
c (double)                  Speed of light in the units used for the simulation.
source_index (int)          Index in the `particles` array for the massive central body.
=========================== ==================================================================

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

**Effect Structure**: *rebx_params_gr_full*

=========================== ==================================================================
Field (C type)              Description
=========================== ==================================================================
c (double)                  Speed of light in the units used for the simulation.
=========================== ==================================================================

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
Python Example          `Migration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Migration.ipynb>`_,
                        `EccAndIncDamping.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/EccAndIncDamping.ipynb>`_.
======================= ===============================================

This updates particles' positions and velocities between timesteps to achieve the desired changes to the osculating orbital elements (exponential growth/decay for a, e, inc, linear progression/regression for Omega/omega.
This nicely isolates changes to particular osculating elements, making it easier to interpret the resulting dynamics.  
One can also adjust the coupling parameter `p` between eccentricity and semimajor axis evolution, as well as whether the damping is done on Jacobi, barycentric or heliocentric elements.
Since this method changes osculating (i.e., two-body) elements, it can give unphysical results in highly perturbed systems.

**Effect Structure**: *rebx_params_modify_orbits_direct*

=========================== ==================================================================
Field (C type)              Description
=========================== ==================================================================
p (double)                  Coupling parameter between eccentricity and semimajor axis evolution
                            (see Deck & Batygin 2015). `p=0` corresponds to no coupling, `p=1` to
                            eccentricity evolution at constant angular momentum.
coordinates (enum)          Type of elements to use for modification (Jacobi, barycentric or heliocentric).
                            See the examples for usage.
=========================== ==================================================================

**Particle Parameters**

=========================== =========== ======================================================
Name (C type)               Required    Description
=========================== =========== ======================================================
tau_a (double)              No          Semimajor axis exponential growth/damping timescale
tau_e (double)              No          Eccentricity exponential growth/damping timescale
tau_inc (double)            No          Inclination axis exponential growth/damping timescale
tau_Omega (double)          No          Period of linear nodal precession/regression
tau_omega (double)          No          Period of linear apsidal precession/regression
=========================== =========== ======================================================

.. _modify_orbits_forces:

modify_orbits_forces
********************

======================= ===============================================
Authors                 D. Tamayo, H. Rein
Implementation Paper    *In progress*
Based on                `Papaloizou & Larwood 2000 <http://labs.adsabs.harvard.edu/adsabs/abs/2000MNRAS.315..823P/>`_.
C Example               :ref:`c_example_modify_orbits`
Python Example          `Migration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Migration.ipynb>`_
                        `EccAndIncDamping.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/EccAndIncDamping.ipynb>`_.
======================= ===============================================

This applies physical forces that orbit-average to give exponential growth/decay of the semimajor axis, eccentricity and inclination.
The eccentricity damping keeps the angular momentum constant (corresponding to `p=1` in modify_orbits_direct), which means that eccentricity damping will induce some semimajor axis evolution.
Additionally, eccentricity/inclination damping will induce pericenter/nodal precession.
Both these effects are physical, and the method is more robust for strongly perturbed systems.

**Effect Structure**: *rebx_params_modify_orbits_forces*

=========================== ==================================================================
Field (C type)              Description
=========================== ==================================================================
coordinates (enum)          Type of elements to use for modification (Jacobi, barycentric or heliocentric).
                            See the examples for usage.
=========================== ==================================================================

**Particle Parameters**

=========================== =========== ======================================================
Name (C type)               Required    Description
=========================== =========== ======================================================
tau_a (double)              No          Semimajor axis exponential growth/damping timescale
tau_e (double)              No          Eccentricity exponential growth/damping timescale
tau_inc (double)            No          Inclination axis exponential growth/damping timescale
tau_Omega (double)          No          Period of linear nodal precession/regression
tau_omega (double)          No          Period of linear apsidal precession/regression
=========================== =========== ======================================================

Radiation Forces
^^^^^^^^^^^^^^^^

.. _radiation_forces:

radiation_forces
****************

======================= ===============================================
Authors                 H. Rein, D. Tamayo
Implementation Paper    *In progress*
Based on                `Burns et al. 1979 <http://labs.adsabs.harvard.edu/adsabs/abs/1979Icar...40....1B/>`_.
C Example               :ref:`c_example_rad_forces_debris_disk`, :ref:`c_example_rad_forces_circumplanetary`.
Python Example          `Radiation_Forces_Debris_Disk.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Radiation_Forces_Debris_Disk.ipynb>`_,
                        `Radiation_Forces_Circumplanetary_Dust.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Radiation_Forces_Circumplanetary_Dust.ipynb>`_.
======================= ===============================================

This applies radiation forces to particles in the simulation.  
It incorporates both radiation pressure and Poynting-Robertson drag.
Only particles whose `beta` parameter is set will feel the radiation.  

**Effect Structure**: *rebx_params_radiation_forces*

=========================== ==================================================================
Field (C type)              Description
=========================== ==================================================================
c (double)                  Speed of light in the units used for the simulation.
source_index (int)          Index in the `particles` array for the radiation source.
=========================== ==================================================================

**Particle Parameters**

=========================== =========== ======================================================
Name (C type)               Required    Description
=========================== =========== ======================================================
beta (double)               No          Ratio of the radiation force to the gravitational force
                                        from the radiation source.
=========================== =========== ======================================================

.. _paramlist:

Table of Particle Parameters
----------------------------


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

