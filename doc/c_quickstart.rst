.. _c_quickstart:

Quickstart Guide (C)
====================

Installation
------------

Navigate to the parent directory that holds the ``rebound`` folder (see below if you want to install in a custom folder).  Then::

    git clone https://github.com/dtamayo/reboundx.git

(install git if you don't have it).  *If you do use* a custom install location for REBOUNDx, you have to additionally set the ``REB_DIR`` environment variable to the path to REBOUND. You might add this to your shell's startup files, e.g. with bash,::
    
    export REB_DIR=/Users/dtamayo/rebound

.. _c_qs:

Quick Start Guide
-----------------

We assume we've already set up a ``reb_simulation`` called ``sim``.  We always begin by adding REBOUNDx to our simulation::
    
    struct rebx_extras* rebx = rebx_init(sim);

We then add the effect we are interested in, which returns a pointer to a rebx_effect to which we can later add parameters.
For example, to add general relativity effects,::

    struct rebx_effect* gr = rebx_add_effect(rebx, "gr");

For details on particular effect implementations, and for a list of available effects, see :ref:`effects`.

We can now set parameters to our particles and effects.  
**To understand which parameters need to be set for particular effects, consult** :ref:`effects`.
For example, general relativity effects require us to set the speed of light in appropriate units.

We always first add a parameter to the effect (or particle), specifying the appropriate type from the rebx_param_type enumeration (see :ref:`c_api` under enums).
This gives us back a pointer that we can then use to set the value::

    double* c = rebx_add_param(gr, "c", REBX_TYPE_DOUBLE);
    *c = 3.e8;

The same goes for a particle, e.g. setting a semimajor axis damping timescale for the modify_orbits_forces effect::

    double* tau_a = rebx_add_param(&sim->particles[1], "tau_a", REBX_TYPE_DOUBLE);
    *tau_a = -1.e4;

Note that rebx_add_param and rebx_get_param (below) always take a pointer to the particle or the effect.

We can check the value of a particle's parameter with the pointer we got back (c and tau_a above), but if particles or effects are removed from the simulation, dereferencing these pointers leads to undefined behavior / segmentation faults.
When accessing parameters, one should therefore retrieve them each time with rebx_get_param::

    double* tau_a = rebx_get_param(&sim->particles[1], "tau_a");

This allows you to check whether the pointer is valid, as rebx_get_param will return NULL if the parameter is not found (**so make sure you check for NULL!**).

You can add as many effects as you'd like in the same simulation.
Once you're done setting up all the effects you want, just run the REBOUND simulation as usual::

    reb_integrate(sim, tmax);

Probably the quickest way to get up and running is to modify an existing example for your effect.
You can find links to the appropriate examples here: :ref:`effects`, along with implementation descriptions and citations.
You can find the example files in the ``reboundx/examples`` folder.

Even if you are using the C version, you might also take a look at the python example links at :ref:`effects`, as the iPython notebooks nicely incorporate text and may therefore have a bit longer discussions about the physical details for each effect.

**Note that REBOUNDx is an all-or-nothing proposition.  Either you use it for all additional effects, or none.  
If you use REBOUNDx and then try to set sim->additional_forces with your own custom routine, you will overwrite the function pointer that REBOUNDx is using under the hood.
We therefore provide functions for adding your own custom effects within REBOUNDx itself.
Follow the :ref:`c_example_custom_ptm` tutorial for how to do this.**  
