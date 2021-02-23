.. _c_quickstart:

Quickstart (C)
==============

Installation
------------

Navigate to the parent directory that holds the ``rebound`` folder (see below if you want to install in a custom folder).  Then in a terminal::

    git clone https://github.com/dtamayo/reboundx.git

(install git if you don't have it).  *If you do use* a custom install location for REBOUNDx, you have to additionally set the ``REB_DIR`` environment variable to the path to REBOUND. You might add this to your shell's startup files, e.g. with bash::
    
    export REB_DIR=/Users/dtamayo/rebound

.. _c_qs:

Quick Start Guide
-----------------

We assume we've already set up a ``reb_simulation`` called ``sim``.  We always begin by attaching REBOUNDx to our simulation::
    
    struct rebx_extras* rebx = rebx_attach(sim);

We then add the effect we are interested in.
There are two types of effects, forces and operators. 
The easiest is to check the documentation for the effect you're interested in on the :ref:`effects` page, which links to a C example that shows you how to add it and set the relevant parameters.
For example, to add post-Newtonian forces::

    struct rebx_force* gr = rebx_load_force(rebx, "gr");

This returns a rebx_force pointer. 
Different effects have their own set of parameters that must be set, listed on :ref:`effects` and explained in the examples. 
For example, general relativity effects require us to set the speed of light (always in the units we're using in the rest of the simulation).

Simple parameters of type double can be set like this::

    rebx_set_param_double(rebx, &gr->ap, "c", 3e4);

Additional parameters (ap) are stored in a linked list for each force, operator and particle, so in addition to passing
the rebx instance, we always pass a pointer to the head of the ap linked list (&gr->ap). We then pass the name of 
the parameter, and the value we want to set it to.

The same goes for a particle. For example, to set a semimajor axis damping timescale tau_a for the modify_orbits_forces effect to particles[1]::

    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -1.e4);

We can check the value of a particle's parameter by getting a pointer to that parameter with::

    double* c = rebx_get_param(rebx, gr->ap, "c");
   
Note that if the parameter name (here "c") is not found, rebx_get_param will return NULL.
So check for NULL before dereferencing it, or you will get undefined behavior/seg faults.

Once you've set up your force, you still have to add it to the simulation. You do that with::

    rebx_add_force(rebx, force);

You can add as many effects as you'd like in the same simulation.
Once you're done setting up all the effects you want, just run the REBOUND simulation as usual::

    reb_integrate(sim, tmax);

Probably the quickest way to get up and running is to modify an existing example for your effect.
You can find links to the appropriate examples here: :ref:`effects`, along with implementation descriptions.
You can find the example files in the ``reboundx/examples`` folder.

Even if you are using the C version, you might also take a look at the python example links at :ref:`effects`, as the iPython notebooks nicely incorporate text and may therefore have a bit longer discussions about the physical details for each effect.
