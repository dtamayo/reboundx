.. _python_quickstart:

Quickstart (Python)
===================

Installation
------------

To quickly get up and running, simply type into a terminal 
(if you installed REBOUND in a virtual environment, activate it first)::

    pip install reboundx

At this point you are done and can skip to the Quick Start Guide below.

For a more complete installation, i.e., if you want any of the following: 

* Source code
* The example files so that you can modify them locally.
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
``reboundx.Extras`` instance, passing it the simulation we want to attach it to:

.. code:: python

    import reboundx
    rebx = reboundx.Extras(sim)

We then add the effect we are interested in.
There are two types of effects, forces and operators. 
The easiest is to go to the :ref:`effects` page, which lists all effects and links to a jupyter
notebook example of how to set it up.
For a deeper discussion, see Tamayo et al. 2019.
Loading a force/operator returns an object of the appropriate type.
For example, let's add some mass loss to the star:

.. code:: python

    mm = rebx.load_operator("modify_mass")
    rebx.add_operator(mm)

Each effect will have different parameters to set, listed on the :ref:`effects` page and the examples.
Forces, operators and particles have a params attribute that works like a dictionary.
For example, let's add an exponential mass loss (i.e., negative) timescale to the star (index 0) of 1000 time units.

.. code:: python

    sim.particles[0].params["tau_mass"] = -1000

You can add as many modifications as you'd like in the same simulation.
Simply add them:

.. code:: python

    gr = rebx.load_force("gr")
    rebx.add_force(gr)
    gr.params['c'] = 1.e4 # set speed of light

The units for the various parameters should always match the units you're using for the rest of the simulation (see the examples).
When you're done setting up the modifications you want, you just run your REBOUND simulation like you normally would:

.. code:: python

    sim.integrate(100.)

Probably the quickest way to get up and running is to modify an existing example for your effect.
You can find links to the appropriate examples here: :ref:`effects`, as well as details of each implementation and citations.
