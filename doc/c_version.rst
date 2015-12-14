.. _c_version:

C Version
=========

Installation
------------

Navigate to the parent directory that holds the ``rebound`` folder (see below if you want to install in a custom folder).  Then::

    git clone https://github.com/dtamayo/reboundx.git

(install git if you don't have it).  *If you do use* a custom install location for REBOUNDx, you have to additionally set the ``REB_DIR`` environment variable to the path to REBOUND (you might add this to your shell's startup files)::
    
    export REB_DIR=/Users/dtamayo/rebound

Quick Start Guide
-----------------

We assume we've already set up a ``reb_simulation`` called ``sim``.  We always begin by adding REBOUNDx to our simulation::
    
    struct rebx_extras rebx = rebx_init(sim);

We then add the effect we are interested in::

    rebx_add_effect(rebx);

where ``effect`` is one of the effects in :ref:`effectList`.  Some effects need parameters to set up, see :ref:`add-effects`.  We then set the particle-specific parameters::

    rebx_set_param(&sim->particles[1], 1.e4);

where ``param`` is one of the parameters in :ref:`paramList`.  Here we set the hypothetical ``param`` parameter for ``particles[1]`` to a value of 1.e4.  Then we run the REBOUND simulation as usual::

    reb_integrate(sim, tmax);

The best way to get started is to use the examples as a starting point and modify them as needed.  You can find at least one example for each REBOUNDx effect in the ``reboundx/examples`` folder (see a listing at :ref:`c_examples`).

Even if you are using the C version, you should also take a look at the iPython examples at :ref:`ipython_examples`, as the iPython notebooks nicely incorporate text and they therefore have a bit longer discussions about the physical details for each effect.

C API
-----

Top Level REBOUNDx Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygengroup:: MainRebxFunctions

.. _add-effects:

Functions For Adding REBOUNDx Effects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :ref:`effectList` for additional details on the implementation of each effect.

.. doxygengroup:: AddEffect

Functions For Getting and Setting Particle Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both getters and setters always take a *pointer* to the particle.  
See :ref:`paramList` for definitions of the various parameters.

.. doxygengroup:: GettersSetters

Convenience Functions for Particular Effects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are functions provided for convenience for calculating various parameters.

.. doxygengroup:: ConvFunc

Internal Functions and Data Structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For documentation of the code that the user does not diretly interface with, refer to the reboundx.h file.
