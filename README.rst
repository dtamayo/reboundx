.. image:: http://img.shields.io/badge/REBOUNDx-v2.8.4-green.svg?style=flat
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

Update: New Version!
====================

Update: New version!

After incoporating a few different types of effects, we have come up with a more general infrastructure for REBOUNDx.  
The new version allows for:

* Adding effects more than once (e.g. if you wanted to turn tides on or off on different bodies individually)
* Calling REBOUND functions within REBOUNDx 
* Passing messages from REBOUNDx to REBOUND in ipython, so you can see what when wrong rather than the kernel simply dying
* Simpler syntax when adding and changing particle and effect parameters
* More robust automatic installation with pip

This should allow for a stable API moving forward as people add new effects.
You should therefore update both REBOUND and REBOUNDx to the latest versions (see Sec. 5.3 of http://rebound.readthedocs.org/en/latest/python_quickstart.html)
Let me know if you have any feedback / issues.

Welcome to REBOUNDx (eXtras)
====================================

REBOUNDx allows you to easily incorporate additional physics into your REBOUND simulations.
All the computationally expensive parts of the code are written in C, so that the code will run much faster than if you define your own custom `additional_forces` functions in Python.

For a list of supported effects, installation instructions, tutorials/examples and documentation, please see http://reboundx.readthedocs.org.

.. image:: https://github.com/dtamayo/dtamayo.github.io/blob/master/pix/reboundx.png

