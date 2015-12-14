.. _python_version:

Python Version
==============

Installation
------------

To quickly get up and running, simply type into a terminal 
(if you installed REBOUND in a virtual environment, activate it first)::

    pip install reboundx

If you don't have REBOUND installed, this command will also install it for you.
At this point you are done and can skip to the Quick Start Guide below.

For a more complete installation, i.e., if you want any of the following: 

* Source code
* Example files that you can modify (you can always view the examples on `Github
  <https://github.com/dtamayo/reboundx/tree/master/ipython_examples/>`_)
* To also use the C version
  
Navigate to the parent directory that holds the ``rebound`` folder and 
(see below if you want to install in a custom folder)::

    git clone https://github.com/dtamayo/reboundx.git

(install git if you don't have it).  
To install the Python version from this repository 
(you'd also do this to install the Python version after modifying any of the C code)::

    cd reboundx
    pip install -e ./

To use a custom install location for REBOUNDx, you have to set the ``REB_DIR`` environment 
variable to the path to REBOUND (you might add this to your shell's startup files)::
    
    export REB_DIR=/Users/dtamayo/rebound

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
``reboundx.Extras`` instance, passing it the simulation we want to
modify:

.. code:: python

    import reboundx
    rebx = reboundx.Extras(sim)

We then add the modification we are interested in. 
We do this with member functions that all follow the same recipe: 
``add_`` + the name of the modification (see :ref:`effectList`). 
For example:

.. code:: python

    rebx.add_modify_orbits_forces()

In general, for each effect there are also particle-specific parameters. 
These are set directly through the particle:

.. code:: python

    sim.particles[1].tau_a = -1000.

Here we set the timescale for semimajor axis decay for ``particles[1]``.
Each modification has its own set of parameters.  You can consult the :ref:`paramlist` 
or the appropriate iPython example for the effect you wish to include in :ref:`ipython_examples`.

You can add as many modifications as you'd like in the same simulation.
Simply add them:

.. code:: python

    rebx.add_gr()

When you're done setting up the modifications you want, you just run
your REBOUND simulation like you normally would:

.. code:: python

    sim.integrate(100.)

Check out the other :ref:`ipython_examples` for more concrete examples of the
various modifications.

Python API
----------

.. autoclass:: reboundx.Extras
    :members:
