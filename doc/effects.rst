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

One can pick and choose which particles have which parameters set.  
For each particle, any unset parameter is ignored.

=========================== ======================================================
Name (C type)               Description
=========================== ======================================================
tau_a (double)              Semimajor axis exponential growth/damping timescale
tau_e (double)              Eccentricity exponential growth/damping timescale
tau_inc (double)            Inclination axis exponential growth/damping timescale
tau_Omega (double)          Period of linear nodal precession/regression
tau_omega (double)          Period of linear apsidal precession/regression
=========================== ======================================================

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

One can pick and choose which particles have which parameters set.  
For each particle, any unset parameter is ignored.

=========================== ======================================================
Name (C type)               Description
=========================== ======================================================
tau_a (double)              Semimajor axis exponential growth/damping timescale
tau_e (double)              Eccentricity exponential growth/damping timescale
tau_inc (double)            Inclination axis exponential growth/damping timescale
tau_Omega (double)          Period of linear nodal precession/regression
tau_omega (double)          Period of linear apsidal precession/regression
=========================== ======================================================

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

Only particles with their ``beta`` parameter set will feel radiation forces.

=========================== ======================================================
Name (C type)               Description
=========================== ======================================================
beta (double)               Ratio of the radiation force to the gravitational force
                            from the radiation source.
=========================== ======================================================

Mass modifications
^^^^^^^^^^^^^^^^^^

.. _modify_mass:

modify_mass
***********

Set particles' ``tau_mass`` parameter to a negative value for mass loss, positive for mass growth.

======================= ===============================================
Authors                 D. Tamayo
Implementation Paper    *In progress*
Based on                None
C Example               :ref:`c_example_modify_mass`
Python Example          `ModifyMass.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/ModifyMass.ipynb>`_.
======================= ===============================================

This adds exponential mass growth/loss to individual particles every timestep.

**Effect Structure**:

*None*

**Particle Parameters**

Only particles with their ``tau_mass`` parameter set will have their masses affected.

=========================== ======================================================
Name (C type)               Description
=========================== ======================================================
tau_mass (double)           e-folding mass loss (<0) or growth (>0) timescale    
=========================== ======================================================

.. _custom:

custom_force, custom_post_timestep_modification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

======================= ===============================================
Authors                 H. Rein, D. Tamayo
Implementation Paper    N/A
Based on                N/A
C Example               :ref:`c_example_custom_ptm`.
Python Example          N/A
======================= ===============================================

This is in case you want to use your own quick-and-dirty force or post-timestep modification function in your ``problem.c`` file (i.e. with the C version).
If you want to use it in Python, you might as well set up a new effect in the REBOUNDx framework (see the add an effect section).
You can also write your own quick-and-dirty Python function, but this will switch between C and Python every timestep and will be slower than implementing the effect in C by a factor of a few.
See `Forces.ipynb <https://github.com/hannorein/rebound/blob/master/ipython_examples/Forces.ipynb>`_ for how to do this.
