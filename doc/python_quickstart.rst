.. _python_quickstart:

Quickstart Guide (Python)
=========================

Installation
------------

To quickly get up and running, simply type into a terminal 
(if you installed REBOUND in a virtual environment, activate it first)::

    pip install reboundx

At this point you are done and can skip to the Quick Start Guide below.

For a more complete installation, i.e., if you want any of the following: 

* Source code
* The example files that so you can modify them locally.
* To also use the C version
 
First follow the installation instructions for the C version in :ref:`c_quickstart`.
Then, to install the Python version from this repository, navigate to the `reboundx` directory and
(you'd also do this to install the Python version after modifying any of the C code)::

    pip install -e ./

.. _python_qs:

Quick Start Guide
-----------------

You always begin by setting up your REBOUND simulation like you normally
would, e.g.

.. code:: python

    import rebound
    sim = rebound.Simulation()
    sim.add(m=1.)
    sim.add(a=1.)

To use reboundx, we first import it, and then create a
``reboundx.Extras`` instance, passing it the simulation we want to modify:

.. code:: python

    import reboundx
    rebx = reboundx.Extras(sim)

We then add the modification we are interested in. 
We do this with member functions that all follow the same naming recipe: 
``add_`` + the name of the modification (see :ref:`effects`).
Each of these functions returns an instance of the parameters class appropriate for the added effect.
For example:

.. code:: python

    params = rebx.add_modify_orbits_forces()

We can then set parameters for the effect as a whole directly, e.g.,

.. code:: python

    params.p = 0.5

This set the coupling parameter between eccentricity and semimajor axis evolution for all particles.
In general, for each effect there are also particle-specific parameters. 
These are set directly through the particle:

.. code:: python

    sim.particles[1].tau_a = -1000.
    print(sim.particles[1].tau_a)

Here we set the timescale for semimajor axis decay for ``particles[1]``.
In general, each effect has its own particular set of parameters, both for the effect as a whole, and for individual particles.

**The main reference point in the documentation is** :ref:`effects` **,which has descriptions for each effect and its parameters, citations, and links to examples.**

You can find descriptions of each effect's adder function and any convenience functions at :ref:`python_api`.

You can add as many modifications as you'd like in the same simulation.
Simply add them:

.. code:: python

    rebx.add_gr()

You can even add the same effect more than once if you want.
When you're done setting up the modifications you want, you just run your REBOUND simulation like you normally would:

.. code:: python

    sim.integrate(100.)

Probably the quickest way to get up and running is to modify an existing example for your effect.
You can find links to the appropriate examples here: :ref:`effects`.
