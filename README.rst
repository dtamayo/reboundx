Welcome to REBOUNDx (eXtras)
====================================

.. image:: https://github.com/dtamayo/dtamayo.github.io/blob/master/pix/reboundx.png

REBOUNDx allows you to easily incorporate additional physics into your REBOUND simulations.
All the computationally expensive parts of the code are written in C, so that the code will run much faster than if you define your own custom `additional_forces` functions in Python.

So far we include:

* General relativity corrections
* Semimajor axis and eccentricity damping, implemented as forces
* Direct modifications to particles' orbital elements after each timestep
* Radiation Forces

Installation
------------

You can call REBOUNDx from whichever language you use for REBOUND (C or Python).

If you want to quickly get up and running with the *Python version*, simply type into a terminal (if you installed REBOUND in a virtual environment, activate it first)::

    pip install reboundx

If you don't have REBOUND installed, this command will also install it for you.

For a more complete installation, i.e., if you want any of the following: 
* Source code
* Example files
* To use the *C version* (or both C and Python)
  
Navigate to the parent directory that holds the `rebound` folder and (see below if you want to install in a custom folder)::

    git clone https://github.com/dtamayo/reboundx.git

(install git if you don't have it).  You can now run C code (see the rebound/examples directory).  To install the Python version from this repository (you'd also do this to install the Python version after modifying any of the C code)::

    cd reboundx
    pip install -e ./
   
If you cloned the repository, and want to use a custom install location for REBOUNDx, you have to set the `REB_DIR` environment variable to the path to REBOUND.  You might add this to your shell's startup files, e.g. .bashrc or .profile::
    
    export REB_DIR=/Users/dtamayo/rebound

Getting Started
---------------

The best way to get started is to use the examples as a starting point and modify them as needed.  Even if you didn't clone the repository, you can still see the examples at https://github.com/dtamayo/reboundx/tree/master/ipython_examples (for the Python examples) and https://github.com/dtamayo/reboundx/tree/master/examples (for the C examples).  If you are using C, you might still look through the corresponding Python example, as the ipython notebooks nicely incorporate text and they therefore have a bit longer discussions about the details for the implementation.

Documentation
-------------

For details, be sure to check out the documentation at http://reboundx.readthedocs.org.
