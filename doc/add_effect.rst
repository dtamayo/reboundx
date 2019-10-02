.. _add_effect:

Adding A New Effect
============================

Whether you're playing around with new physics in your integrations or want to contribute a new effect to REBOUNDx, the process is easy.
Our hope is that as people use the package and work on new problems, they will contribute their new effects so others can use (and cite!) their implementations.
Contributing the effects you implement (and may want to use in the future) to REBOUNDx has the added benefit of ensuring that they stay up to date as REBOUND and REBOUNDx expand.

*Do I Have To Write It In C?*

Forces and operators are called every timestep, and the overhead of REBOUND calling a Python function each timestep makes it a factor of a few slower than if the effect was written in C.
Therefore all effects in REBOUNDx are written in C.
Often you might want to quickly try something out in Python, and the Custom Effects ipynb example shows you how to do that.
Once you want to do something more serious, or contribute to REBOUNDx, starting from the skeleton of one of the other effects and following this guide is straightforward, even if you don't normally code in C.

*Writing Forces*

There are two kinds of effects, forces and operators.
We will go through the process of adding a force in detail, and comment on the very similar process for operators and the relevant differences at the end.
Adding a force boils down to writing a C function that will evaluate and add the relevant accelerations to all the particles.
Then behind the scenes, REBOUNDx then takes care of when your acceleration function gets called depending on the REBOUND integrator the user has chosen.

Here we take writing the ``radiation_forces.c`` force as an example, adding features in steps.

Note: Different implementations for the same effect should be added as separate effects.
For example, there are separate files for ``gr``, ``gr_potential``, and ``gr_full``.

Writing the C Code
------------------ 

effect.c
^^^^^^^^^^^^^^^^^^^^^

First copy an existing ``.c`` file in ``reboundx/src``.
You should copy existing ones from other effects, so that you have the license and right code structure to work from.
Typically you'll be able to reuse many parts of the code.
Here we erase everything below the comments and write it from scratch for clarity.

In our case, our function reads

.. code-block:: c

    rebx_radiation_forces(struct reb_simulation* const sim, struct rebx_force* const radiation_forces, struct reb_particle* const particles, const int N){
        double c = 3.e8;
        double beta = 0.1;
        struct reb_particle source = particles[0];

        const double mu = sim->G*source.m;
        for (int i=0;i<N;i++){
            
            if(i == source_index) continue;
            
            const struct reb_particle p = particles[i];
            const double dx = p.x - source.x; 
            const double dy = p.y - source.y;
            const double dz = p.z - source.z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz); // distance to star
            
            const double dvx = p.vx - source.vx;
            const double dvy = p.vy - source.vy;
            const double dvz = p.vz - source.vz;
            const double rdot = (dx*dvx + dy*dvy + dz*dvz)/dr; // radial velocity
            const double a_rad = *beta*mu/(dr*dr);

            // Equation (5) of Burns, Lamy & Soter (1979)

            particles[i].ax += a_rad*((1.-rdot/c)*dx/dr - dvx/c);
            particles[i].ay += a_rad*((1.-rdot/c)*dy/dr - dvy/c);
            particles[i].az += a_rad*((1.-rdot/c)*dz/dr - dvz/c);
        }
    }

This gives a recipe for accessing the effect parameters, checking for NULL to avoid segmentation faults and warn the user when needed parameters are not set.
It's also important to use ``_N_real`` for the number of particles in the simulation, since ``sim->N`` includes any variational particles that have been added.
The idea is to search for a particle that has a ``radiation_source`` flag set, and then call an internal function that does the force calculation given the source's index in the ``particles`` array.
This is the safe way for the user to use the effect, since it is robust against particles moving around in the ``particles`` array.
If particle has the ``radiation_source`` parameter set, it defaults to the 0 index particle being the source, which makes sense for simple cases.

Contact me if you need to add support for parameters with different types than those in the rebx_param_type enumeration (see the Enums section of :ref:`c_api`).

core.c and core.h
^^^^^^^^^^^^^^^^^

You need to add your new effect in the rebx_add function in reboundx/src/core.c.
It should be self-explanatory to mirror what other effects are doing, but make sure you set *either* effect->force or effect->ptm depending on whether your effect is a force or post_timestep_modification.
Additionally, if you are implementing a force, and your force depends on the particle velocities, you need to include a sim->force_is_velocity_dependent = 1 line like, e.g., gr_full.
You also need to add your function prototype at the bottom of reboundx/src/core.h.

Python Code
^^^^^^^^^^^

Particle Parameters
^^^^^^^^^^^^^^^^^^^
Particularly if you are used to the Python side of REBOUND/REBOUNDx, you should read :ref:`c_quickstart` to see how to access parameters, and consult effects that are already written.

What parameters would you want to set for individual particles?

In our case, we can imagine different particles having different sizes, which affects the acceleration they feel from the radiation.
One approach would be to set all the relevant physical properties of each particle, e.g., bulk density, physical radius, radiation pressure coefficient etc.
We instead chose to use the single dimensionless quantity ``beta`` needed for the force calculation, and add convenience functions for calculating ``beta`` from physical parameters (see below).
This means a bit more setup for the user, but a cleaner and more efficient implementation that doesn't have to look up several parameters for each particle in performance-critical code.
I don't know the right answer, but it is probably worth spending a bit of time thinking about these design decisions.

The only constraint is that you cannot use particle parameter names that are in use by other effects, unless they certainly refer to the same property, e.g., ``bulk_density`` would probably be fine even if another effect used it. 
You can check the master list of all parameter names currently used in REBOUNDx, with links to what they mean, here: :ref:`params`.

Effect Parameters
^^^^^^^^^^^^^^^^^

What parameters would you want to set for the effect as a whole?

REBOUND allows the user to set whatever set of units they want to use for their simulations.  
So for any parameter with dimensions that is required by your effect, you will need the user to pass the appropriate value.
In our case, we need to know the speed of light.

Additionally, we choose to allow the user to specify which particle in the simulation is the radiation source.


Example/Test Case
^^^^^^^^^^^^^^^^^

All effects have a corresponding example (typically adapted from code to test the implementation) that others can work from.

Navigate to the ``reboundx/examples`` folder, and copy the ``modify_orbits`` folder to another folder named after your effect.

We now also want to update all the ``Makefiles`` and setup scripts to include your new effect.
If you navigate to ``reboundx/scripts`` and type ``python add_new_effect.py``, the script will automatically detect the new effect file and make all the required changes.

Go back to ``reboundx/examples/youreffect/`` and modify ``problem.c`` file as you like.
You can then run your program in your example folder, typing ``make`` (you should  ``make clean`` first if you make changes to the code in reboundx/src), and then ``./rebound``.
All examples use a standard Makefile that compiles and links all the required libraries, so you shouldn't have to edit it.  

If you get an error about OpenGL or GLUT, just google `install openGL glut libraries <your OS here>` for instructions, or open your ``Makefile`` and set OPENGL=0 (it's easier to debug if you can see what's going on though!)
See Sec. 2.4 of `OpenGL Keyboard Commands <http://rebound.readthedocs.org/en/latest/c_quickstart.html>`_ for a list of the visualization keyboard commands.

iPython Example
^^^^^^^^^^^^^^^

If you don't use iPython notebooks, you should try them!
I use them for all my (research) dynamics simulations.
All the Python examples in REBOUND and REBOUNDx also use them.
iPython is now part of the Jupyter project, and you can find installation instructions `here <http://jupyter.readthedocs.org/en/latest/install.html>`_.

I think most people using REBOUND/REBOUNDx use the Python implementation, so if you're up for it, add an iPython notebook in ``reboundx/ipython_examples/``.
You might copy ``EccAndIncDamping.ipynb`` and edit that as a starter.

Add Your Effect to the Main Documentation Page!
-----------------------------------------------

You add the documentation for your effect directly within your ``effect.c`` file.
It will then automatically get built into the :ref:`effects` page.
Easiest is if you copy-paste from another effect source file.

At the top of the cmoment block, you should edit the file, brief and author lines.
The rest of the documentation goes Below the dollar signs.
In our case, 
 
.. code-block:: rst

     * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     *
     * $Radiation Forces$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
     *
     * ======================= ===============================================
     * Authors                 H. Rein, D. Tamayo
     * Implementation Paper    *In progress*
     * Based on                `Burns et al. 1979 <http://labs.adsabs.harvard.edu/adsabs/abs/1979Icar...40....1B/>`_.
     * C Example               :ref:`c_example_rad_forces_debris_disk`, :ref:`c_example_rad_forces_circumplanetary`.
     * Python Example          `Radiation_Forces_Debris_Disk.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Radiation_Forces_Debris_Disk.ipynb>`_,
     *                         `Radiation_Forces_Circumplanetary_Dust.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Radiation_Forces_Circumplanetary_Dust.ipynb>`_.
     * ======================= ===============================================
     * 
     * This applies radiation forces to particles in the simulation.  
     * It incorporates both radiation pressure and Poynting-Robertson drag.
     * Only particles whose `beta` parameter is set will feel the radiation.  
     * 
     * **Effect Parameters**
     * 
     * ============================ =========== ==================================================================
     * Field (C type)               Required    Description
     * ============================ =========== ==================================================================
     * c (double)                   Yes         Speed of light in the units used for the simulation.
     * ============================ =========== ==================================================================
     *
     * **Particle Parameters**
     *
     * If no particles have radiation_source set, effect will assume the particle at index 0 in the particles array is the source.
     *
     * ============================ =========== ==================================================================
     * Field (C type)               Required    Description
     * ============================ =========== ==================================================================
     * radiation_source (int)       No          Flag identifying the particle as the source of radiation.
     * ============================ =========== ==================================================================
     * 
     */

We first add the group that our effect belongs to, between dollar signs, $Radiation Forces$.
This keeps different implementations of, e.g., general relativity corrections in the same place.
If you want to make a new category for your effect, edit :ref:`effect_headers` (/reboundx/doc/effect_headers.rst).
You can optionally add a description general to all implementations in the category following the format in the file, which will show up in :ref:`effects`.

Then fill in the table:
``Authors`` says who wrote the code.
``Implementation paper`` is the paper that you'd like to be cited by people using your implementation.
``Based on`` is the paper that the equations you used come from.

``C Example`` is a link to the C Example you wrote.
All C examples in the ``reboundx/examples`` directory are automatically built into the documentation, and have cross-reference targets of the form ``c_example_foldername``, where foldername is the name of your example folder in ``reboundx/examples``. 

For the ``Python Example`` line, edit the link from another documentation entry with the name of your ipython notebook filename (in both the title and bracketed URL).

Underneath your table, provide a description that will inform users when it's appropriate to apply your effect (and when it's not!).

Finally, if your effect requires the user to set (possibly optionally) particular effect or particle parameters, we create tables for them too. 

You can check how everything looks by navigating to ``reboundx/doc`` and typing ``make clean``, then ``make html``.
Then navigate to ``reboundx/doc/_build/html`` and open ``index.html`` in your browser.
The main effects page (with the tables) is on the left: REBx Effects & Parameters.
The automatically included documentation will be under API Documentation (Python) and API Documentation (C).

.. _pullrequest:

Putting together a Pull Request
-------------------------------

If you'd rather e-mail me your code, I'm happy to incorporate it, but if you'd like for github to show your account as a contributor to the project, send me a pull request! 

If you have never used git, it's very useful for backups, rewinding errors, and collaboration.
You can make an account at `http://github.com <http://github.com>`_.
Follow the instructions under `Time to Submit Your First PR` `here <http://www.thinkful.com/learn/github-pull-request-tutorial/Expect-a-Thorough-Review#Time-to-Submit-Your-First-PR>`_ up until "Tadaa!" to fork the REBOUNDx repository and make your own local branch.

Now you can modify the code as described below, and can incrementally commit changes.
As a starting point, you can check out `this guide <https://www.atlassian.com/git/tutorials/saving-changes>`_.

After working through this document and making all the changes, you can then send me a pull request by following the rest of the instructions in the pull request tutorial above.

reboundx.h (Optional)
^^^^^^^^^^^^^^^^^^^^^

If you want to provide any convenience functions for the user, add the prototypes at the bottom under ``Convenience functions for various effects``.
Include some mention of your effect (in short form) in the function name, and follow the format for other functions to have the documentation automatically built into reboundx.readthedocs.org.
In our case

.. code-block:: c

    /**
     * @brief Calculates beta, the ratio between the radiation pressure force and the gravitational force from the star.
     * @param G Gravitational constant.
     * @param c Speed of light.
     * @param source_mass Mass of the source body.
     * @param source_luminosity Luminosity of radiation source.
     * @param radius Particle physical radius.
     * @param density density of particle.
     * @param Q_pr Radiation pressure coefficient (Burns et al. 1979).
     * @return Beta parameter (double). 
     */
    double rebx_rad_calc_beta(const double G, const double c, const double source_mass, const double source_luminosity, const double radius, const double density, const double Q_pr);

Python Code
-----------

With the REBOUNDx version, your effect will automatically work from Python.
You only have to add a couple lines of code if you added a convenience function for the user, or if you defined new structures for your particular effect.
I'm happy to help with the latter.

First navigate to ``reboundx/`` and type ``pip install -e .``.
This will install the updated libreboundx extension so you can call it from Python.
You'll have to run the same command any time you edit the C code (you don't need to after changing the Python code--if using an ipython notebook, just restart the kernel after making changes to the Python code).

Now open ``reboundx/reboundx/extras.py``.

Following our example:

.. code-block:: python

    def rad_calc_beta(self, G, c, source_mass, source_luminosity, radius, density, Q_pr):
        """
        Calculates a particle's beta parameter (the ratio of the radiation force to the gravitational force).
        All values must be passed in the same units as used for the simulation as a whole (e.g., AU, Msun, yr/2pi).

        :param G: Gravitational constant
        :param c: Speed of light
        :param source_mass: Mass of radiation source
        :param source_luminosity: Luminosity of radiation source
        :param radius: grain's physical radius
        :param density: particle bulk density
        :param Q_pr: radiation pressure coefficient
        :type G: float
        :type c: float
        :type source_mass: float
        :type source_luminosity: float
        :type radius: float
        :type density: float
        :type Q_pr: float
        :rtype: float
        """
        clibreboundx.rebx_rad_calc_beta.restype = c_double
        return clibreboundx.rebx_rad_calc_beta(c_double(G), c_double(c), c_double(source_mass), c_double(source_luminosity), c_double(radius), c_double(density), c_double(Q_pr))

The documentation is self-explanatory (just follow same format), and as above will get automatically built into the online documentation.
In the code, the first line tells ``ctypes`` what to expect the C function to return (here a ``double``).
In the last line, we again cast everything to ``ctypes`` types, and for any parameters the C function expects as a pointer, we use ``byref()``.
See the ctypes documentation for details: https://docs.python.org/3/library/ctypes.html or contact me for help.

An operator is a function that, given a timestep, updates all the particle states (i.e., the positions, velocities, masses, or any other parameters you want to add to particles).
One reason to use these would be when evolving things other than the particles in the simulation, e.g., a disk model.
This is the idea of splitting methods like the Wisdom-Holman map, and can yield very efficient integration schemes.
If you're interested in doing this, you should read the REBOUNDx paper.
Good examples to go off are the ``modify_mass`` and ``modify_orbits_direct`` effects.

