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
    sim.add(m=1., hash="sun")
    sim.add(a=1., hash="earth")

To use reboundx, we first import it, and then create a
``reboundx.Extras`` instance, passing it the simulation we want to modify:

.. code:: python

    import reboundx
    rebx = reboundx.Extras(sim)

We then add the modification we are interested in (for a listing see :ref"`effects`).
Each of these functions returns an instance of the rebx_effect class.
For example:

.. code:: python

    effect = rebx.add("modify_orbits_direct")

Effects and particles have a params attribute that works like a dictionary.
We set parameters for the effect with

.. code:: python

    effect.params["p"] = 0.5

This sets the coupling parameter between eccentricity and semimajor axis evolution for all particles.
In general, for each effect there are also particle-specific parameters. 
These are set similarly:

.. code:: python
    sim.particles["earth"].params["tau_a"] = -1000.

Here we set the timescale for semimajor axis decay for ``earth``.
In general, each effect has its own particular set of parameters, both for the effect as a whole, and for individual particles.

**The user is responsible for checking what parameters need to be set for a particular effect, which can be found at** :ref:`effects`.

You can add as many modifications as you'd like in the same simulation.
Simply add them:

.. code:: python

    rebx.add("gr")

When you're done setting up the modifications you want, you just run your REBOUND simulation like you normally would:

.. code:: python

    sim.integrate(100.)

Probably the quickest way to get up and running is to modify an existing example for your effect.
You can find links to the appropriate examples here: :ref:`effects`, as well as details of each implementation and citations.
