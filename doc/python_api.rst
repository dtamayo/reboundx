.. _python_api:

API Documentation (Python)
==========================

This documents the adder functions for each REBOUNDx effect, including default behaviors, and any convenience functions for particular effects.
**The main reference point in the documentation is** :ref:`effects` **,which has descriptions for each effect and its parameters, citations, and links to examples.**

.. automodule:: reboundx
    :members:

Convenience Functions
---------------------

There are a number of convenience functions (for calculating the potential associated with a particular conservative effect, or converting between sets of parameters), which are listed in the "Convenience Functions for Various Effects" section of :ref:`c_api`.  
For example, rebx_gr_potential_hamiltonian() should be called as rebx.gr_potential_hamiltonian() from Python.
In all cases, the pointer to the simulation is what you get from rebound.Simulation(), and pointers to rebx_effects are what you get back from the corresponding rebx.add("effect").
These convenience functions are typically shown in the corresponding Python example for a particular effect. 

