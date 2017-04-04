.. image:: http://img.shields.io/badge/REBOUNDx-v2.17.2-green.svg?style=flat
    :target: http://reboundx.readthedocs.org
.. image:: https://badge.fury.io/py/reboundx.svg
    :target: https://badge.fury.io/py/reboundx
.. image:: https://travis-ci.org/dtamayo/reboundx.svg?branch=master
    :target: https://travis-ci.org/dtamayo/reboundx
.. image:: https://coveralls.io/repos/dtamayo/reboundx/badge.svg?branch=master&service=github 
    :target: https://coveralls.io/github/dtamayo/reboundx?branch=master
.. image:: http://img.shields.io/badge/license-GPL-green.svg?style=flat 
    :target: https://github.com/dtamayo/reboundx/blob/master/LICENSE
.. image:: https://readthedocs.org/projects/pip/badge/?version=latest
    :target: http://reboundx.readthedocs.org/
.. image:: https://img.shields.io/badge/launch-binder-ff69b4.svg?style=flat
    :target: http://mybinder.org/repo/dtamayo/reboundx

 Welcome to REBOUNDx (eXtras)
====================================

REBOUNDx allows you to easily incorporate additional physics into your REBOUND simulations.
All the computationally expensive parts of the code are written in C, so that the code will run much faster than if you define your own custom `additional_forces` functions in Python.

For a list of supported effects, installation instructions, tutorials/examples and documentation, please see http://reboundx.readthedocs.org.

.. image:: https://github.com/dtamayo/dtamayo.github.io/blob/master/pix/reboundx.png

Changelog
=========
  - 2.17.0 Added support for applying forces as leapfrog operators (for velocity-dependent forces)
  - 2.16.0 Added support for binaries and simulationarchive
  - 2.14.0 Added track_min_distance
  - 2.13.0 Fixed collision for long integrations and added support for arrays
  - 2.12.0 Added central force effect
  - 2.11.0 Added tidal precession
  - 2.10.0 Streamlined code for adding new effects (see below)

Update: REBOUNDx 2.16.0
=======================

We have added support for writing and loading REBOUNDx binary save files, making it possible to use REBOUND SimulationArchives for machine-independent, bit-by-bit reproducible simulations with additional effects. We also improved the support for array parameters, and updated the documentation. We added information in the binary files of the githash used to generate it, making it easier to match up versions to reproduce simulations.

Update: REBOUNDx 2.13.0
=======================

We have added support for array parameters, and have settled on what we think should be a general enough API going forward.
We also fixed a memory collision issue for long integrations, so you should update to the new version when possible.

Update: REBOUNDx 2.10.0
=======================

After getting feedback from users and some experimentation with different types of effects one might wish to add to simulations, we have implemented a more streamlined scheme for adding new effects.  
In this new version, code that you add in C will be automatically callable from the Python version without any of the coding overhead that was required before.
For details, see http://reboundx.readthedocs.io/en/latest/add_effect.html

On the Python side We have also developed a custom dictionary for accessing/setting effect/particle parameters (before we were using a quick and dirty solution of hijacking default getters and setters).
This has changed some of the syntax, and previous code using REBOUNDx will need (small) modifications to work with version 2.10.0.
The changes should be clear from inspection of the various examples in the documentation (see links in http://reboundx.readthedocs.io/en/latest/effects.html).
Some of these changes led to feature additions in REBOUND, so you should also update REBOUND to at least version 2.19.2 (see Sec. 5.3 of http://rebound.readthedocs.org/en/latest/python_quickstart.html).
Let me know if you have any feedback / issues.


