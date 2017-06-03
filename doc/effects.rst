.. _effects:

REBx Effects & Parameters
=========================

Below are descriptions for each of the effects included in REBOUNDx.
Different implementations for the same effect are grouped together.
All effects follow the same recipes for usage, see the Python quick-start guide (:ref:`python_qs`) or C quick-start guide (:ref:`c_qs`).
Probably the quickest way to get up and running is to edit one of the linked examples for the effect you're interested in.

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

**Effect Parameters**

If p is not set, it defaults to 1.  If coordinates not set, defaults to using Jacobi coordinates.

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
p (double)                   No          Coupling parameter between eccentricity and semimajor axis evolution
                                         (see Deck & Batygin 2015). `p=0` corresponds to no coupling, `p=1` to
                                         eccentricity evolution at constant angular momentum.
coordinates (enum)           No          Type of elements to use for modification (Jacobi, barycentric or particle).
                                         See the examples for usage.
============================ =========== ==================================================================

**Particle Parameters**

One can pick and choose which particles have which parameters set.  
For each particle, any unset parameter is ignored.

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
tau_a (double)               No          Semimajor axis exponential growth/damping timescale
tau_e (double)               No          Eccentricity exponential growth/damping timescale
tau_inc (double)             No          Inclination axis exponential growth/damping timescale
tau_Omega (double)           No          Period of linear nodal precession/regression
tau_omega (double)           No          Period of linear apsidal precession/regression
============================ =========== ==================================================================


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

**Effect Parameters**

If coordinates not, defaults to using Jacobi coordinates.

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
coordinates (enum)           No          Type of elements to use for modification (Jacobi, barycentric or particle).
                                         See the examples for usage.
============================ =========== ==================================================================

**Particle Parameters**

One can pick and choose which particles have which parameters set.  
For each particle, any unset parameter is ignored.

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
tau_a (double)               No          Semimajor axis exponential growth/damping timescale
tau_e (double)               No          Eccentricity exponential growth/damping timescale
tau_inc (double)             No          Inclination axis exponential growth/damping timescale
============================ =========== ==================================================================


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
Adding this effect to several bodies is NOT equivalent to using gr_full.

**Effect Parameters**

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
c (double)                   Yes         Speed of light in the units used for the simulation.
============================ =========== ==================================================================

**Particle Parameters**

If no particles have gr_source set, effect will assume the particle at index 0 in the particles array is the source.

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
gr_source (int)              No          Flag identifying the particle as the source of perturbations.
============================ =========== ==================================================================


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

This algorithm incorporates the first-order post-newtonian effects from all bodies in the system, and is necessary for multiple massive bodies like stellar binaries.

**Effect Parameters**

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
c (double)                   Yes         Speed of light in the units used for the simulation.
============================ =========== ==================================================================

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

**Effect Parameters**

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
c (double)                   Yes         Speed of light in the units used for the simulation.
============================ =========== ==================================================================

**Particle Parameters**

If no particles have gr_source set, effect will assume the particle at index 0 in the particles array is the source.

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
gr_source (int)              No          Flag identifying the particle as the source of perturbations.
============================ =========== ==================================================================


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

**Effect Parameters**

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
c (double)                   Yes         Speed of light in the units used for the simulation.
============================ =========== ==================================================================

**Particle Parameters**

If no particles have radiation_source set, effect will assume the particle at index 0 in the particles array is the source.

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
radiation_source (int)       No          Flag identifying the particle as the source of radiation.
beta (float)                 Yes         Ratio of radiation pressure force to gravitational force. Particles without beta set feel no radiation forces.
============================ =========== ==================================================================


Mass Modifications
^^^^^^^^^^^^^^^^^^

.. _modify_mass:

modify_mass
***********

======================= ===============================================
Authors                 D. Tamayo
Implementation Paper    *In progress*
Based on                None
C Example               :ref:`c_example_modify_mass`
Python Example          `ModifyMass.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/ModifyMass.ipynb>`_.
======================= ===============================================

This adds exponential mass growth/loss to individual particles every timestep.
Set particles' ``tau_mass`` parameter to a negative value for mass loss, positive for mass growth.

**Effect Parameters**

*None*

**Particle Parameters**

Only particles with their ``tau_mass`` parameter set will have their masses affected.

============================ =========== =======================================================
Name (C type)                Required    Description
============================ =========== =======================================================
tau_mass (double)            Yes         e-folding mass loss (<0) or growth (>0) timescale    
============================ =========== =======================================================


Tides
^^^^^^^^^^^^^^^^^^

.. _tides_precession:

tides_precession
****************

======================= ===============================================
Authors                 D. Tamayo
Implementation Paper    *In progress*
Based on                `Hut 1981 <https://ui.adsabs.harvard.edu/#abs/1981A&A....99..126H/abstract>`_.
C Example               :ref:`c_example_tides_precession`.
Python Example          `TidesPrecession.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesPrecession.ipynb>`_.
======================= ===============================================

This adds precession from the tidal interactions between the particles in the simulation and the central body, both from tides raised on the primary and on the other bodies.
In all cases, we need to set masses for all the particles that will feel these tidal forces. After that, we can choose to include tides raised on the primary, on the "planets", or both, by setting the respective bodies' R_tides (physical radius) and k1 (apsidal motion constant, half the tidal Love number).
You can specify the primary with a "primary" flag.
If not set, the primary will default to the particle at the 0 index in the particles array.

**Effect Parameters**

None

**Particle Parameters**

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
R_tides (float)              Yes         Physical radius (required for contribution from tides raised on the body).
k1 (float)                   Yes         Apsidal motion constant (half the tidal Love number k2).
primary (int)                No          Set to 1 to specify the primary.  Defaults to treating particles[0] as primary if not set.
============================ =========== ==================================================================


.. _tides_synchronous_ecc_damping:

tides_synchronous_ecc_damping
*****************************

======================= ===============================================
Authors                 D. Tamayo
Implementation Paper    *In progress*
Based on                `Hut 1981 <https://ui.adsabs.harvard.edu/#abs/1981A&A....99..126H/abstract>`_.
C Example               :ref:``.
Python Example          ` <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesPrecession.ipynb>`_.
======================= ===============================================

This adds precession from the tidal interactions between the particles in the simulation and the central body, both from tides raised on the primary and on the other bodies.
In all cases, we need to set masses for all the particles that will feel these tidal forces. After that, we can choose to include tides raised on the primary, on the "planets", or both, by setting the respective bodies' R_tides (physical radius) and k1 (apsidal motion constant, half the tidal Love number).
You can specify the primary with a "primary" flag.
If not set, the primary will default to the particle at the 0 index in the particles array.

**Effect Parameters**

None

**Particle Parameters**

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
R_tides (float)              Yes         Physical radius (required for contribution from tides raised on the body).
k1 (float)                   Yes         Apsidal motion constant (half the tidal Love number k2).
primary (int)                No          Set to 1 to specify the primary.  Defaults to treating particles[0] as primary if not set.
============================ =========== ==================================================================


Central Force
^^^^^^^^^^^^^^^^^^

.. _central_force:

central_force
*************

======================= ===============================================
Authors                 D. Tamayo
Implementation Paper    *In progress*
Based on                None
C Example               :ref:`c_example_central_force`
Python Example          `CentralForce.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/CentralForce.ipynb>`_.
======================= ===============================================

Adds a general central acceleration of the form a=Aradial*r^gammaradial, outward along the direction from a central particle to the body.
Effect is turned on by adding Aradial and gammaradial parameters to a particle, which will act as the central body for the effect,
and will act on all other particles.

**Effect Parameters**

None

**Particle Parameters**

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
Aradial (double)             Yes         Normalization for central acceleration.
gammaradial (double)         Yes         Power index for central acceleration.
============================ =========== ==================================================================


Gravity Fields
^^^^^^^^^^^^^^^^^^

.. _gravitational_harmonics:

gravitational_harmonics
***********************

======================= ===============================================
Authors                 D. Tamayo
Implementation Paper    *In progress*
Based on                None
C Example               :ref:`c_example_J2`
Python Example          `J2.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/J2.ipynb>`_.
======================= ===============================================


**Effect Parameters**

None

**Particle Parameters**

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
J2 (double)                  No          J2 coefficient
J4 (double)                  No          J4 coefficient
R_eq (double)                Yes         Equatorial radius of nonspherical body used for calculating Jn harmonics
============================ =========== ==================================================================


Miscellaneous Utilities
^^^^^^^^^^^^^^^^^^^^^^^

.. _track_min_distance:

track_min_distance
******************

======================= ===============================================
Authors                 D. Tamayo
Implementation Paper    *In progress*
Based on                None
C Example               :ref:`c_example_track_min_distance`
Python Example          `TrackMinDistance.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TrackMinDistance.ipynb>`_.
======================= ===============================================

For a given particle, this keeps track of that particle's minimum distance from another body in the simulation.  User
should add parameters to the particular particle whose distance should be tracked.

**Effect Parameters**

*None*

**Particle Parameters**

Only particles with their ``min_distance`` parameter set initially will track their minimum distance. The effect will
update this parameter when the particle gets closer than the value of ``min_distance``, so the user has to set it
initially.  By default distance is measured from sim->particles[0], but you can specify a different particle by setting
the ``min_distance_from`` parameter to the hash of the target particle.

================================ =========== =======================================================
Name (C type)                    Required    Description
================================ =========== =======================================================
min_distance (double)            Yes         Particle's miminimum distance.
min_distance_from (uint32)       No          Hash for particle from which to measure distance
min_distance_orbit (reb_orbit)   No          Parameter to store orbital elements at moment corresponding to min_distance (heliocentric)
================================ =========== =======================================================


