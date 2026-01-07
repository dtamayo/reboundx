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

*None*

**Particle Parameters**

*None*


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


