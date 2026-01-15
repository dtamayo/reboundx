.. _collisions:

Collisions
==========

Below are each of the functions for resolving collisions implemented in REBOUNDx.

.. _fragmenting_collisions:

fragmenting_collisions
**********************

======================= ===============================================
Authors                 H. Tajer, H. Rein, T. Lu
Implementation Paper    `In prep.`_.
Based on                `Leinhardt and Stewart 2012 <https://iopscience.iop.org/article/10.1088/0004-637X/745/1/79>`_, `Chambers 2013 <https://www.sciencedirect.com/science/article/pii/S0019103513000754?via%3Dihub>`_.
C Example               :ref:`c_example_fragmenting_collisions_embryo_disk`
Python Example          `FragmentingCollisions.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/FragmentingCollisions.ipynb>`_.
======================= ===============================================

This module is based on the collision prescription of Leinhardt & Stewart (2012), from now on called LS2012. 
After a collision between two bodies is detected, the outcome will depend on the velocity of the collision, mass of the bodies, and the impact angle.
Possible outcomes are:

1. "merging": particles will merge 

2. "fragmentation": erosion of the target or accretion of projectile material into target (fragmentation)

3. "elastic bounce": particles will leave each other intact

For more details, refer to the implementation paper (in prep).

**Effect Parameters**

======================================== =========== ==============================================================================
Field (C type)                           Required    Description
======================================== =========== ==============================================================================
fc_min_frag_mass (double)                 Yes         Minimum fragment mass allowed in the simulation. 
fc_separation_distance_scale (double)     No          Ratio of distance between the COM of newly added fragment and COM of target and projectile over sum of radii of target and projectile
fc_rho1 (double)                          No          Density rho_1 used in computing the spherical radius of combined mass of target and projectile if they had density rho_1. Used in computing the disruption criteria for objects with different bulk densities and mass ratios. For more information, refer to Leinhardt and Stewart 2012.
fc_cstar (double)                         No          C* is a dimentionless parameter used as a measure of the dissipation of energy within the target. For more information refer to Leinhardt and Stewart 2012.
fc_particle_list_file (string)            No          Name of the output file. Example is: "family_tree.csv". If you wish to produce the output file, you need to set this parameter. 
fc_id_max (int)                           No          Maximum particle ID assigned so far. This parameter gets updated everytime a new particle ID is assigned.
======================================== =========== ==============================================================================

**Particle Parameters**

======================================== =========== ==============================================================================
Field (C type)                           Required    Description
======================================== =========== ==============================================================================
fc_id (int)                              No          Unique particle ID. Everytime a collision happens, new particle IDs are assigned to all the bodies after the collision, and previous IDs are discarded.
======================================== =========== ==============================================================================


.. _merging_collisions:

merging_collisions
******************

======================= ===============================================
Authors                 H. Rein 
Based on                None
C Example               :ref:`c_example_merging_collisions`
======================= ===============================================

This is a simple example implementation of a REBOUNDx collision module.
The outcome is similar to the built-in REBOUND function reb_collision_resolve_merge.

**Effect Parameters**

*None*

**Particle Parameters**

*None*


