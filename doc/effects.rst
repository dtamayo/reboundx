.. _effects:

Implemented Effects
===================

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

.. _modify_orbits_forces:

modify_orbits_forces
********************

======================= ===============================================
Authors                 D. Tamayo, H. Rein
Implementation Paper    `Kostov et al., 2016 <https://ui.adsabs.harvard.edu/abs/2016ApJ...832..183K/abstract>`_.
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


.. _modify_orbits_direct:

modify_orbits_direct
********************

======================= ===============================================
Authors                 D. Tamayo
Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_. 
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

If p is not set, it defaults to 0.  If coordinates not set, defaults to using Jacobi coordinates.

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


.. _exponential_migration:

exponential_migration
*********************

======================= ===============================================
Author                   Mohamad Ali-Dib
Implementation Paper    `Ali-Dib et al., 2021 AJ <https://arxiv.org/abs/2104.04271>`_.
Based on                `Hahn & Malhotra 2005 <https://ui.adsabs.harvard.edu/abs/2005AJ....130.2392H/abstract>`_.
C Example               :ref:`c_example_exponential_migration`
Python Example          `ExponentialMigration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/ExponentialMigration.ipynb>`_.
======================= ===============================================

Continuous velocity kicks leading to exponential change in the object's semimajor axis. 
One of the standard prescriptions often used in Neptune migration & Kuiper Belt formation models.
Does not directly affect the eccentricity or inclination of the object.

**Particle Parameters**

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
em_tau_a (double)              Yes          Semimajor axis exponential growth/damping timescale
em_aini (double)               Yes          Object's initial semimajor axis
em_afin (double)               Yes          Object's final semimajor axis
============================ =========== ==================================================================


General Relativity
^^^^^^^^^^^^^^^^^^

.. _gr_potential:

gr_potential
************

======================= ===============================================
Authors                 H. Rein, D. Tamayo
Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
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


.. _gr:

gr
**

======================= ===============================================
Authors                 P. Shi, D. Tamayo, H. Rein
Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
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
Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
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

Radiation Forces
^^^^^^^^^^^^^^^^

.. _radiation_forces:

radiation_forces
****************

======================= ===============================================
Authors                 H. Rein, D. Tamayo
Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
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
Implementation Paper    `Kostov et al., 2016 <https://ui.adsabs.harvard.edu/abs/2016ApJ...832..183K/abstract>`_.
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

.. _tides_constant_time_lag:

tides_constant_time_lag
***********************

======================= ===============================================
Authors                 Stanley A. Baronett, D. Tamayo, Noah Ferich
Implementation Paper    Baronett et al., in prep.
Based on                `Hut 1981 <https://ui.adsabs.harvard.edu/#abs/1981A&A....99..126H/abstract>`_, `Bolmont et al., 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...583A.116B/abstract>`_.
C Example               :ref:`c_example_tides_constant_time_lag`.
Python Example          `TidesConstantTimeLag.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesConstantTimeLag.ipynb>`_.
======================= ===============================================

This adds constant time lag tidal interactions between orbiting bodies in the simulation and the primary, both from tides raised on the primary and on the other bodies.
In all cases, we need to set masses for all the particles that will feel these tidal forces. After that, we can choose to include tides raised on the primary, on the "planets", or both, by setting the respective bodies' physical radius particles[i].r, k1 (apsidal motion constant, half the tidal Love number), constant time lag tau, and rotation rate Omega. See Hut (1981) and Bolmont et al. 2015 above.

If tau is not set, it will default to zero and yield the conservative piece of the tidal potential.

**Effect Parameters**

None

**Particle Parameters**

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
particles[i].r (float)       Yes         Physical radius (required for contribution from tides raised on the body).
tctl_k1 (float)                   Yes         Apsidal motion constant (half the tidal Love number k2).
tctl_tau (float)                  No          Constant time lag. If not set will default to 0 and give conservative tidal potential
Omega (float)                No          Rotation rate. If not set will default to 0
============================ =========== ==================================================================


Central Force
^^^^^^^^^^^^^^^^^^

.. _central_force:

central_force
*************

======================= ===============================================
Authors                 D. Tamayo
Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
Based on                None
C Example               :ref:`c_example_central_force`
Python Example          `CentralForce.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/CentralForce.ipynb>`_.
======================= ===============================================

Adds a general central acceleration of the form a=Acentral*r^gammacentral, outward along the direction from a central particle to the body.
Effect is turned on by adding Acentral and gammacentral parameters to a particle, which will act as the central body for the effect,
and will act on all other particles.

**Effect Parameters**

None

**Particle Parameters**

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
Acentral (double)             Yes         Normalization for central acceleration.
gammacentral (double)         Yes         Power index for central acceleration.
============================ =========== ==================================================================


Gravity Fields
^^^^^^^^^^^^^^^^^^

.. _gravitational_harmonics:

gravitational_harmonics
***********************

======================= ===============================================
Authors                 D. Tamayo
Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_. 
Based on                None
C Example               :ref:`c_example_J2`
Python Example          `J2.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/J2.ipynb>`_.
======================= ===============================================

Adds azimuthally symmetric gravitational harmonics (J2, J4) to bodies in the simulation. Current implementation assumes everything is planar, i.e. spin pole of body aligned with z axis of simulation.

**Effect Parameters**

None

**Particle Parameters**

============================ =========== ==================================================================
Field (C type)               Required    Description
============================ =========== ==================================================================
J2 (double)                  No          J2 coefficient
J4 (double)                  No          J4 coefficient
R_eq (double)                No          Equatorial radius of nonspherical body used for calculating Jn harmonics
============================ =========== ==================================================================


Integration Steppers
^^^^^^^^^^^^^^^^^^^^

These are wrapper functions to taking steps with several of REBOUND's integrators in order to build custom splitting schemes.

.. _steppers:

steppers
********

======================= ===============================================
Authors                 D. Tamayo, H. Rein
Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
Based on                `Rein and Liu, 2012 <https://ui.adsabs.harvard.edu/abs/2012A%26A...537A.128R/abstract>`_.
C Example               None
Python Example          `CustomSplittingIntegrationSchemes.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/CustomSplittingIntegrationSchemes.ipynb>`_.
======================= ===============================================

These are wrapper functions to taking steps with several of REBOUND's integrators in order to build custom splitting schemes.

**Effect Parameters**

None

**Particle Parameters**

None


Parameter Interpolation
^^^^^^^^^^^^^^^^^^^^^^^

This isn't an effect that's loaded like the others, but an object that facilitates machine-independent interpolation of parameters that can be shared by both the C and Python versions. See the examples below for how to use them.

.. _interpolation:

interpolation
*************

======================= ===============================================
Authors                 S.A. Baronett, D. Tamayo, N. Ferich
Implementation Paper    Baronett et al., in prep.
Based on                `Press et al., 1992 <https://ui.adsabs.harvard.edu/abs/1992nrca.book.....P/abstract>`_. 
C Example               :ref:`c_example_parameter_interpolation`
Python Example          `ParameterInterpolation.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/ParameterInterpolation.ipynb>`_.
======================= ===============================================

**Effect Parameters**

Not applicable. See examples.

**Particle Parameters**

Not applicable. See examples.

Miscellaneous Utilities
^^^^^^^^^^^^^^^^^^^^^^^

.. _track_min_distance:

track_min_distance
******************

======================= ===============================================
Authors                 D. Tamayo
Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
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
initially.  By default, distance is measured from sim->particles[0], but you can specify a different particle by setting
the ``min_distance_from`` parameter to the hash of the target particle.

================================ =========== =======================================================
Name (C type)                    Required    Description
================================ =========== =======================================================
min_distance (double)            Yes         Particle's mininimum distance.
min_distance_from (uint32)       No          Hash for particle from which to measure distance
min_distance_orbit (reb_orbit)   No          Parameter to store orbital elements at moment corresponding to min_distance (heliocentric)
================================ =========== =======================================================


