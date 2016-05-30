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

We then add the effect we are interested in, which returns a pointer to a parameter struct that you can later modify::

    struct rebx_params_effect* params = rebx_add_effect(rebx);

where ``effect`` is one of the effects in :ref:`effects`.
Some effects adders require you to pass parameter values, see :ref:`c_api`.
We can now set effect parameters directly, e.g.,::

    params->c = 3.e8;

We then set particle-specific parameters with::

    rebx_set_param_double(&sim->particles[1], "tau_a", 1.e4);

Note that
    * There is a different setter function for each type, named ``rebx_set_param_type`` (replace type, as above).
    * We always pass the address to a particle in a simulation.
    * We pass the name of the parameter we want to change.

We can then check the value of a particle's parameter with::

    rebx_get_param_double(&sim->particles[1], "tau_a");

In general, each effect has its own particular set of parameters, both for the effect as a whole, and for individual particles.

**The main reference point in the documentation is** :ref:`effects` **,which has descriptions for each effect and its parameters, citations, and links to examples.**

You can find descriptions of each effect's adder function and any convenience functions at :ref:`c_api`.


You can add as many modifications as you'd like in the same simulation (even the same effect more than once).
Once you're done setting up all the modifications you want, just run the REBOUND simulation as usual::

    reb_integrate(sim, tmax);

Probably the quickest way to get up and running is to modify an existing example for your effect.
You can find links to the appropriate examples here: :ref:`effects`.
You can find the example files in the ``reboundx/examples`` folder.

Even if you are using the C version, you might also take a look at the python example links at :ref:`effects`, as the iPython notebooks nicely incorporate text and they therefore have a bit longer discussions about the physical details for each effect.

**Note that REBOUNDx is an all-or-nothing proposition.  Either you use it for all additional effects, or none.  
If you use REBOUNDx and then try to set sim->additional_forces with your own custom routine, you will overwrite the function pointer that REBOUNDx is using under the hood.
We therefore provide functions for adding your own custom effects within REBOUNDx itself.
Follow the :ref:`c_example_custom_ptm` tutorial for how to do this.**  
