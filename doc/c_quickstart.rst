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
We can now set effect parameters.  
For example, general relativity effects require us to set the speed of light in appropriate units,  e.g.,::

    rebx_set_effect_param_double(gr, "c", 3.e8);

Notice that we specify in the function name that we are adding to an `effect` and we specify the variable type.  
We can use an analogous function to set particle-specific parameters, e.g.::

    rebx_set_particle_param_double(&sim->particles[1], "tau_a", 1.e4);

Note that
    * Getter and setter functions always take a pointer to the particle or the effect.
    * We pass the name of the parameter we want to change.

We can then check the value of a particle's parameter with (replace `particle` with `effect` for checking an effect's parameters)::

    double* tau_a = rebx_get_particle_param_double(&sim->particles[1], "tau_a");

**Important:** If the parameter is not found (e.g., if we asked for "tau_e"), the tau_a pointer we get back will be set to NULL, so if we try to dereference it we will get a segmentation fault.
It is therefore prudent to always check the return value for NULL before dereferencing.::

In general a particular effect will require you to set particle-specific parameters in particles, and general effect parameters in the effect.
**The user is responsible for checking what parameters need to be set for a particular effect, which can be found at** :ref:`effects`.

You can add as many modifications as you'd like in the same simulation.
Once you're done setting up all the modifications you want, just run the REBOUND simulation as usual::

    reb_integrate(sim, tmax);

Probably the quickest way to get up and running is to modify an existing example for your effect.
You can find links to the appropriate examples here: :ref:`effects`, along with implementation descriptions and citations.
You can find the example files in the ``reboundx/examples`` folder.

Even if you are using the C version, you might also take a look at the python example links at :ref:`effects`, as the iPython notebooks nicely incorporate text and may therefore have a bit longer discussions about the physical details for each effect.

**Note that REBOUNDx is an all-or-nothing proposition.  Either you use it for all additional effects, or none.  
If you use REBOUNDx and then try to set sim->additional_forces with your own custom routine, you will overwrite the function pointer that REBOUNDx is using under the hood.
We therefore provide functions for adding your own custom effects within REBOUNDx itself.
Follow the :ref:`c_example_custom_ptm` tutorial for how to do this.**  
