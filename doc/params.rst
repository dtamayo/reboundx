.. _params:

Table of Particle Parameters
============================

This is a master list of all the particle parameter names used by all effects in REBOUNDx.
When adding a new effect, one must avoid name clashes with any other effect that may reasonably be added at the same time as your new one, unless they would definitely refer to the same property.
For example, it would be fine for `bulk_density` to be used for two different effects, but you couldn't use `tau_a` for your a particle's time in the AGB phase in your new effect, since `tau_a` is already taken for semimajor axis evolution.

=========================== ========================================= 
Parameter name              Effects                                  
=========================== ========================================= 
tau_a                       :ref:`modify_orbits_direct`, :ref:`modify_orbits_forces`
tau_e                       :ref:`modify_orbits_direct`, :ref:`modify_orbits_forces`
tau_inc                     :ref:`modify_orbits_direct`, :ref:`modify_orbits_forces`
tau_Omega                   :ref:`modify_orbits_direct`, :ref:`modify_orbits_forces`
tau_omega                   :ref:`modify_orbits_direct`, :ref:`modify_orbits_forces`
beta                        :ref:`radiation_forces`
tau_mass                    :ref:`modify_mass`
=========================== ========================================= 


