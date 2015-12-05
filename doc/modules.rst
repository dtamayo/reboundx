REBOUNDx Effects and Parameters
===============================

.. _effectList:

List of REBOUNDx Effects
------------------------

Each of these can be added to a REBOUND simulation with the corresponding function in :ref:`add-effects`.  Each of these effects has at least one corresponding examples in both ``reboundx/examples`` and ``reboundx/ipython_examples``.  The latter can also be viewed online at https://github.com/dtamayo/reboundx/tree/master/ipython_examples.  Make sure to look at the corresponding iPython example, which has further details on the physical model for the particular effect.

To see what particle parameters should be set for each effect, see :ref:`paramlist`.

=======================  ============================================ 
Effect                   Description
=======================  ============================================ 
modify_orbits_direct      Update orbital elements directly between timesteps
modify_orbits_forces      Modify orbital elements with additional forces
gr                        Post-Newtonian corrections
REB_GRAVITY_TREE          Oct tree, Barnes & Hut 1986, O(N log(N))
REB_GRAVITY_OPENCL        (upgrade to REBOUND 2.0 still in progress) Direct summation, O(N^2), but accelerated using the OpenCL framework.
REB_GRAVITY_FFT           (upgrade to REBOUND 2.0 still in progress) Two dimensional gravity solver using FFTW, works in a periodic box and the shearing sheet. 
=======================  ============================================ 

.. _paramlist:

List of Particle Parameters
---------------------------

=======================  ============================================ 
Module name               Description
=======================  ============================================ 
REB_COLLISION_NONE        No collision detection, default
REB_COLLISION_DIRECT      Direct nearest neighbour search, O(N^2)
REB_COLLISION_TREE        Oct tree, O(N log(N))
REB_COLLISION_SWEPP       (upgrade to REBOUND 2.0 still in progress) Plane sweep algorithm, ideal for low dimensional  problems, O(N) or O(N^1.5) depending on geometry 
=======================  ============================================ 


Boundary conditions
-------------------

=======================  ============================================ 
Module name               Description
=======================  ============================================ 
REB_BOUNDARY_NONE         Dummy. Particles are not affected by boundary conditions, default
REB_BOUNDARY_OPEN         Particles are removed from the simulation if they leaves the box.
REB_BOUNDARY_PERIODIC     Periodic boundary conditions. Particles are reinserted on the other side if they cross the box boundaries. You can use an arbitrary number of ghost-boxes with this module.
REB_BOUNDARY_SHEAR        Shear periodic boundary conditions. Similar to periodic boundary conditions, but ghost-boxes are moving with constant speed, set by the shear.
=======================  ============================================ 
 

Integrators
-----------

=======================  ============================================ 
Module name               Description
=======================  ============================================ 
REB_INTEGRATOR_IAS15      IAS15 stands for Integrator with Adaptive Step-size control, 15th order. It is a vey high order, non-symplectic integrator which can handle arbitrary (velocity dependent) forces and is in most cases accurate down to machine precision. IAS15 can integrate variational equations. Rein & Spiegel 2015, Everhart 1985, default
REB_INTEGRATOR_WHFAST     WHFast is the integrator described in Rein & Tamayo 2015, it's a second order symplectic Wisdom Holman integrator with 11th order symplectic correctors. It is extremely fast and accurate, uses Gauss f and g functions to solve the Kepler motion and can integrate variational equations.
REB_INTEGRATOR_EULER      Euler scheme, first order
REB_INTEGRATOR_LEAPFROG   Leap frog, second order, symplectic
REB_INTEGRATOR_WH         SWIFT-style Wisdom-Holman Mapping, mixed variable symplectic integrator for the Kepler potential, second order, note that  `integrator_whfast.c` almost always offers better characteristics, Wisdom & Holman 1991, Kinoshita et al 1991
REB_INTEGRATOR_SEI        Symplectic Epicycle Integrator (SEI), mixed variable symplectic integrator for the shearing sheet, second order, Rein & Tremaine 2011
REB_INTEGRATOR_HYBRID     An experimental hybrid symplectic integrator that uses WHFast for long term integrations but switches over to IAS15 for close encounters.
=======================  ============================================ 


