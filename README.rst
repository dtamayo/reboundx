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

Installation
------------

You'll have to have git installed.  In a terminal, navigate to the location where you want to put the reboundx folder.  Then::

    git clone https://github.com/dtamayo/reboundx.git

REBOUNDx also needs to know where REBOUND is installed on your machine.  You tell it by setting the REB_DIR environment variable.  One good way to do this is to edit your .profile::

    cd ~
    pico .profile

and add the line::

    export REB_DIR=/Users/dtamayo/rebound

where you would replace the path to your own installation.  Now close your terminal, open a new one and you're all set!
    
Examples
--------

You can find working examples in reboundx/examples.  The gr example adds post-newtonian corrections. The modify_orbits example shows how to add semimajor axis / eccentricity damping to simulations, either through direct modifications of the orbital elements, or through forces.
