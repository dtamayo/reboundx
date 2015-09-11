REBOUNDx (eXtras) 
=================
A library of additional effects/forces for REBOUND simulations
--------------------------------------------------------------

.. image:: https://github.com/dtamayo/dtamayo.github.io/blob/master/pix/reboundx.png

FEATURES
--------

REBOUNDx allows you to easily add typically used modifications to your REBOUND simulations.  So far we include:

* Post-newtonian general relativity corrections
* Semimajor axis and eccentricity damping, implemented as forces
* Direct modifications to particles' orbital elements after each timestep

Python Version Installation
---------------------------

We recommend that you install REBOUND and REBOUNDx within a virtual environment (though this is not necessary).  See http://rebound.readthedocs.org/en/latest/python_quickstart.html.  If you already have REBOUND installed within a virtual environment you don't have to do this step.

Installation is as easy as typing
 
    pip install reboundx

on the command line.  If you don't have REBOUND installed, this command will also install it for you.

You can find examples for adding various effects at https://github.com/dtamayo/reboundx/tree/master/ipython_examples.  You can get these ipython notebooks to play around with them by cloning the directory from github on the command line:

    git clone https://github.com/dtamayo/reboundx.git

(you may have to install git first).  This will create a `reboundx` folder in the current directory.  After you `cd` into `reboundx`, you can install the development version with

    pip install -e ./

(even if you're using a conda environment).  This gives you all the source code and examples, and allows you modify the C code and compile the changes (just type the line above each time you make changes).

C Version Installation
----------------------

In a terminal, navigate to the location where you want to put the reboundx folder.  Then::

    git clone https://github.com/dtamayo/reboundx.git

(you may have to install git first).  REBOUNDx also needs to know where REBOUND is installed on your machine.  You tell it by setting the REB_DIR environment variable.  One good way to do this is to edit your .profile::

    cd ~
    pico .profile

and add the line::

    export REB_DIR=/Users/dtamayo/rebound

where you would replace the path to your own installation.  Now close your terminal, open a new one and you're all set!
    
You can find working examples in reboundx/examples.  The gr example adds post-newtonian corrections. The modify_orbits example shows how to add semimajor axis / eccentricity damping to simulations, either through direct modifications of the orbital elements, or through forces.  It's worth checking out the python examples, since they have a bit more detail about each of the modifications:  https://github.com/dtamayo/reboundx/tree/master/ipython_examples.
