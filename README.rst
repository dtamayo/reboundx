.. image:: http://img.shields.io/badge/REBOUNDx-v3.0.5-green.svg?style=flat
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
.. image:: http://img.shields.io/badge/arXiv-1908.05634-green.svg?style=flat 
    :target: http://arxiv.org/abs/1908.05634
.. image:: https://img.shields.io/badge/launch-binder-ff69b4.svg?style=flat
    :target: http://mybinder.org/repo/dtamayo/reboundx

New Version and Paper!
======================

We've made some big improvements to REBOUNDx to make it more robust and easily extendable, so definitely upgrade to this new version (3.0) if you haven't already.
We've also improved the binary format to better interface with the REBOUND SimulationArchive for the sharing and analysis of machine independent results.

Definitely also check out our paper where we give an overview of the library, analyze the effects of dissipative forces on symplectic integrators, and give some recommendations:`Tamayo, Rein, Shi and Hernandez 2019 <http://arxiv.org/abs/1908.05634>`_

Welcome to REBOUNDx (eXtras)
============================

REBOUNDx allows you to easily incorporate additional physics into your REBOUND simulations.
All the computationally expensive parts of the code are written in C, so that the code will run much faster than if you define your own custom `additional_forces` functions in Python.

For a list of supported effects, installation instructions, tutorials/examples and documentation, please see http://reboundx.readthedocs.io.

.. image:: https://github.com/dtamayo/dtamayo.github.io/blob/master/pix/reboundx.png

Changelog
=========

Version 3.0.5
-------------

* Generalized tides_precession effect to general constant time lag model (Hut 1981)
