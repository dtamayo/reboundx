.. _misc:

Miscellaneous Utilities
=======================

Below are various miscellaneous utilities implemented in REBOUNDx.

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


.. _interpolation:

interpolation
*************

======================= ===============================================
Authors                 S.A. Baronett, D. Tamayo, N. Ferich
Implementation Paper    `Baronett et al., 2022 <https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.6001B/abstract>`_.
Based on                `Press et al., 1992 <https://ui.adsabs.harvard.edu/abs/1992nrca.book.....P/abstract>`_. 
C Example               :ref:`c_example_parameter_interpolation`
Python Example          `ParameterInterpolation.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/ParameterInterpolation.ipynb>`_.
======================= ===============================================

This isn't an effect that's loaded like the others, but an object that facilitates machine-independent interpolation of parameters that can be shared by both the C and Python versions. See the examples for how to use them.

**Effect Parameters**

Not applicable. See examples.

**Particle Parameters**

Not applicable. See examples.

