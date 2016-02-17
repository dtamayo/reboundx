.. _add_effect:

Adding A New Effect
============================

REBOUNDx aims to be modular, so that when adding your own new effect, you don't have to get into all the low-level code details (in ``core.c``).
You just have to write your effect following the REBOUNDx layout, and rely on REBOUNDx to appropriately store, access and call it at the right time.
Here we will code the ``radiation_forces`` effect as an example.

Note: Different implementations for the same effect should be added as separate effects.
For example, there are separate files for ``gr``, ``gr_potential``, and ``gr_full``.
You would then follow these steps for each of these (really not that bad--just copy-paste).

Designing Your Effect
---------------------

REBOUNDx is set up so that you can not only add parameters to individual particles, but also have parameters pertaining to the effect as a whole.  

Particle Parameters
^^^^^^^^^^^^^^^^^^^

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
This adds flexibility for non-standard cases, and makes it possible for the user to add the effect more than once if they for example want to turn on radiation pressure from multiple stars.

After some discussion, we decided that the best way to specify a particle in REBOUND/REBOUNDx is through its index in the ``particles`` array, putting the responsibility on the user to update this effect parameter as needed when particles are added or removed from the simulation.New effects should follow this convention (e-mail me if you don't think this would work for your effect).

.. _pullrequest:

Putting together a Pull Request
-------------------------------

If you'd rather e-mail me your code, I'm happy to incorporate it, but if you'd like for github to show your account as a contributor to the project, send me a pull request! (you'll of course be prominently credited in the documentation no matter how you contribute, see below!).

If you have never used git, it's very useful for backups, rewinding errors, and collaboration.
You can make an account at `http://github.com <http://github.com>`_.
Follow the instructions under `Time to Submit Your First PR` `here <http://www.thinkful.com/learn/github-pull-request-tutorial/Expect-a-Thorough-Review#Time-to-Submit-Your-First-PR>`_ up until "Tadaa!" to fork the REBOUNDx repository and make your own local branch.

Now you can modify the code as described below, and can incrementally commit changes.
As a starting point, you can check out `this guide <https://www.atlassian.com/git/tutorials/saving-changes>`_.

After working through this document and making all the changes, you can then send me a pull request by following the rest of the instructions in the pull request tutorial above.

Writing the C Code
------------------ 

Good documentation is important for new users to actually use your effect. 
REBOUNDx tries to make this process easier by having you document as you write the code, and automating its inclusion into the online documentation.

reboundx.h
^^^^^^^^^^

The place to start is ``reboundx/src/reboundx.h``, which defines the REBOUNDx API (all the functions/data structures the user will use).

First, under ``Parameter structures for each effect``, you should define your effect structure, following the naming convention ``rebx_params_effect``.
If you don't want an effect parameter structure, skip this step.  
In our case

.. code-block:: c

    struct rebx_params_radiation_forces {
         int source_index;                   ///< Index of particle in particles array that is the source of the radiation.
         double c;                           ///< Speed of light in units appropriate for sim->G and initial conditions.
    }

Make sure to use the comment code as above, and it will automatically get built into the documentation!
Then under ``Functions for adding effects``, add the prototype for your effect adder function, which should follow the naming convention ``rebx_add_effect``, should take any effect parameters you want to take from the user, and return a pointer to the structure you defined above (or void with no parameter structure).
In our case

.. code-block:: c

    /**
    * @brief Adds radiation forces to the simulation (i.e., radiation pressure and Poynting-Robertson drag).
    * @param rebx pointer to the rebx_extras instance
    * @param source_index Index in the particles array of the body that is the radiation source.
    * @param c Speed of light.
    */
    struct rebx_params_radiation_forces* rebx_add_radiation_forces(struct rebx_extras* rebx, int source_index, double c);

Use a comment block following the format above for automatic inclusion into the documentation.  
Keep the brief field to a very basic description of the effect (we'll add the main description below).

Finally, if you want to provide any convenience functions for the user, add the prototypes under ``Convenience functions for various effects`` toward the bottom.
Include some mention of your effect (in very short form!) in the function name.
In our case

.. code-block:: c

    /**
     * @brief Calculates beta, the ratio between the radiation pressure force and the gravitational force from the star.
     * @param rebx pointer to the rebx_extras instance.
     * @param params parameters structure returned when adding effect.
     * @param particle_radius radius of grain.
     * @param density density of particle.
     * @param Q_pr Radiation pressure coefficient (Burns et al. 1979).
     * @param L Luminosity of radiation source.
     */
    double rebx_rad_calc_beta(struct rebx_extras* rebx, struct rebx_params_radiation_forces* params, double particle_radius, double density, double Q_pr, double L);

effect.c and effect.h
^^^^^^^^^^^^^^^^^^^^^

Now we add two new files for your effect in ``reboundx/src``, ``effect.c`` and ``effect.h``.
You should copy existing ones from other effects, so that you have the license and right code structure to work from.
At this point you should ask yourself whether your effect is an additional force or a post timestep modification (i.e., something to do between REBOUND timesteps)?
If you're adding a force, you might copy ``radiation_forces``.
If you're adding a post timestep modification, you might copy ``modify_orbits_direct``.

In our case (``radiation_forces``), we have an additional force, but for example mass loss would be a post timestep modification. 

effect.h
^^^^^^^^

In ``effect.h``, you only have to modify the file, brief, and author fields at the top, the include guards (ifndef, define lines) and substitute the name of your effect in the function name.  
Everything else should be kept exactly the same.

effect.c
^^^^^^^^

In ``effect.c``, we first copy paste the file, brief and author lines from ``effect.h``, and change ``#include effect.h`` from the effect you copied to your new one.
Now we write the effect adder function.
In our case

.. code-block:: c

    struct rebx_params_radiation_forces* rebx_add_radiation_forces(struct rebx_extras* rebx, int source_index, double c){
        struct rebx_params_radiation_forces* params = malloc(sizeof(*params));
        params->c = c;
        params->source_index = source_index;
        int force_is_velocity_dependent = 1;
        rebx_add_force(rebx, params, "radiation_forces", rebx_radiation_forces, force_is_velocity_dependent);
        return params;
    }

First we allocate memory for our parameters structure (just replace ``radiation_forces`` with your own effect name).
Then we initialize the params fields we created for our effect structure (if we made one) with what was passed by the user.
Alternatively, if you think the parameters would rarely be changed, you could set them to a default value, and have the user change the values afterward manually (see e.g., modify_orbits_direct.c).
Then if your force is velocity dependent, set ``force_is_velocity_dependent`` to 1, otherwise to 0.
Finally, leave the ``rebx_add_force`` call the same, just replace ``radiation_forces`` in the two function parameters with your own effect name.

If you're adding a post timestep modification, you don't have to specify ``force_is_velocity_dependent`` (cf. modify_orbits_direct.c).

Finally, if you don't have an effect structure (cf. the ``modify_mass`` effect), you should replace ``params`` with ``NULL`` in the call to ``rebx_add_force`` or ``rebx_add_post_timestep_modification`` (cf. the ``modify_mass`` effect).

Now you have to write the main routine for your effect.
A force would update particles' accelerations, while a post timestep modification would update particles' masses, positions and/or velocities.
You might look at different effect implementations for examples of how to access parameters.
In our case, the top of our function looks like

.. code-block:: c

    void rebx_radiation_forces(struct reb_simulation* const sim, struct rebx_effect* const effect){
        const struct rebx_params_radiation_forces* const params = effect->paramsPtr;
        const double c = params->c;
        const int source_index = params->source_index;
        struct reb_particle* const particles = sim->particles;
        const struct reb_particle source = particles[source_index];
        const double mu = sim->G*source.m;

        const int _N_real = sim->N - sim->N_var;
    #pragma omp parallel for
        for (int i=0;i<_N_real;i++){
            if(i == source_index) continue;
            ... 

This gives a recipe for accessing the effect parameters (effect->paramsPtr is a void pointer so you just have to change the const struct line to rebx_params_effect*).
It's also important to use ``_N_real`` for the number of particles in the simulation, since ``sim->N`` includes any variational particles that have been added.

To access particles' individual parameters, we use

.. code-block:: c
    
    const double beta = rebx_get_param_double(&particles[i], "beta");
    if(isnan(beta)) continue; // only particles with beta set feel radiation forces

One nice feature is that for particle-specific parameters, you (probably) don't have to write any extra code to use them.
REBOUNDx uses a hash function to change the name of your particle parameter (here ``"beta"``) to an integer code, which it uses to set and retrieve particle parameters (stored as a linked list starting from the `ap` void pointer (ap = additional parameters) in rebound.Particle).
For all getting and setting, we always  address (&) to a particle in the simulation.
There also needs to be a separate ``rebx_get_param`` and ``rebx_set_param`` for every different variable type (in this case ``"beta"`` is a double).
You can see which getters and setters are currently implemented here: :ref:`getters`.

The getter for each variable type will return a different default value for cases where the parameter is not set for a particular particle (see :ref:`getters`).
For doubles, it's ``nan``, so when looping over the particles, we check for ``nan`` and skip them if beta is not set.

Example/Test Case
^^^^^^^^^^^^^^^^^

At this point, you're done with the C code, though you might consider testing it!
You can kill two birds with one stone by using your test case as an example that others can work from.

Navigate to the ``reboundx/examples`` folder, and copy the ``modify_orbits`` folder to another folder named after your effect.

We now also want to update all the ``Makefiles`` and setup scripts to include your new effect.
If you navigate to ``reboundx/script`` and type ``python add_new_effect.py``, the script will automatically make all the required changes.

Go back to ``reboundx/examples/youreffect/`` and modify ``problem.c`` file as you like.
You can then run your program by navigating to your example folder, typing ``make`` (you may have to do ``make clean`` and then ``make``), and then ``./rebound``.
All examples use a standard Makefile that compiles and links all the required libraries, so you shouldn't have to edit it.  

If you get an error about OpenGL or GLUT, just google `install openGL glut libraries <your OS here>` for instructions, or open your ``Makefile`` and set OPENGL=0 (it's easier to debug if you can see what's going on though!)
See Sec. 2.4 of `OpenGL Keyboard Commands <http://rebound.readthedocs.org/en/latest/c_quickstart.html>`_ for a list of the visualization keyboard commands.

Writing the Python Code
-----------------------

It's now trivial to make your code callable from Python (even if you don't know Python!).
First navigate to ``reboundx/`` and type ``pip install -e .``.
This will install the updated libreboundx extension so you can call it from Python.
You'll have to run the same command any time you edit the C code (you don't need to after changing the Python code--if using an ipython notebook, just restart the kernel after making changes to the Python code).

Now open ``reboundx/reboundx/extras.py``.

Adder Method
^^^^^^^^^^^^

Under `Functions for adding REBOUNDx effects` you have to add your own effect adder.
Copy paste the ``add_gr_potential`` method and change the name to your effect.
In our case

.. code-block:: python

    def add_radiation_forces(self, source_index=0, c=C_DEFAULT):
        """
        You must pass c (the speed of light) in whatever units you choose if you don't use default units of AU, (yr/2pi) and Msun.
        
        :param source_index: Index in the particles array of the body that is the source of the radiation.
        :param c: Speed of light in appropriate units.
        :type source_index: int
        :type c: float
        :rtype: rebx_params_radiation_forces
        """
        clibreboundx.rebx_add_radiation_forces.restype = POINTER(rebx_params_radiation_forces)
        return clibreboundx.rebx_add_radiation_forces(byref(self), c_int(source_index), c_double(c)).contents

One nice feature of Python is that you can make it optional for the user to pass parameters to the adder method by setting default values (here 0 for source_index and C_DEFAULT for c).

Again, we add comments here in a format that allows them to be automatically incorporated into the online documentation.
The block above ``:param ...`` shows up as a description.
``rtype`` (the return type) should always be your parameters structure, and then you document each parameter with a description (``:param param_name:``) and Python type (``:type param_name:``).

The last two lines call the C library.
In the first, you just have to change the effect name to your own.
In the second, again change the effect name, and then pass the parameters your C adder function needs.
You have to cast all passed parameters to a ``ctypes`` type.
The documentation is `here <https://docs.python.org/2/library/ctypes.html>`_, but you can probably get away with copying what's in other methods (you can also look at the bottom of ``extras.py`` for some more complicated ctypes types).
Contact me for help if needed.

If your effect doesn't have a parameters structure, it's even simpler (see the ``add_modify_mass`` method).

Structure Definition
^^^^^^^^^^^^^^^^^^^^

Still in ``extras.py``, (skip this step if you don't have an effect parameters structure) under `Effect parameter class definitions` you have to define your parameters structure.  
Copy paste an existing definition, and again it's probably enough to figure out the ctypes types from other places in the code.
Make sure that the ``_fields_`` match up exactly with what's in your C structure, and that fields appear in the same order.  
In our case

.. code-block:: python

    class rebx_params_radiation_forces(Structure):
        _fields_ = [("source_index", c_int),
                    ("c", c_double)]

Convenience Functions
^^^^^^^^^^^^^^^^^^^^^

If you created any convenience functions in C, add them under `Convenience Functions` (again in ``extras.py``).
One example:

.. code-block:: python

    def rad_calc_beta(self, params, particle_radius, density, Q_pr, L):
        """
        Calculates a particle's beta parameter (the ratio of the radiation force to the gravitational force).
        All values must be passed in the same units as used for the simulation as a whole (e.g., AU, Msun, yr/2pi).

        :param params: parameters instance returned by add_radiation_forces.
        :param particle_radius: grain's physical radius
        :param density: particle bulk density
        :param Q_pr: radiation pressure coefficient
        :param L: Radiation source's luminosity
        :type params: rebx_params_radiation_forces
        :type particle_radius: float
        :type density: float
        :type Q_pr: float
        :type L: float
        :rtype: float
        """
        clibreboundx.rebx_rad_calc_beta.restype = c_double
        return clibreboundx.rebx_rad_calc_beta(byref(self), byref(params), c_double(particle_radius), c_double(density), c_double(Q_pr), c_double(L))

The documentation works as above.
In the code, the first line tells ``ctypes`` what to expect the C function to return (here a ``double``).
In the last line, we again cast everying to ``ctypes`` types, and for any parameters the C function expects as a pointer, we use ``byref()``.

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

You're done with all the code!
Now you want people to use your effect and cite your breakthrough paper!
The main documentation page that summarizes all the effects in REBOUNDx and provides useful links is ``reboundx/doc/effects.rst``.

The document is divided first into groups.  
For example, ``gr``, ``gr_full`` and ``gr_potential`` are all different implementations for general relativity corrections, and are thus lumped under the `General Relativity` heading.  
If your effect fell into an existing category, you would put it there.
In our case, there are no other radiation forces implementations, so we start a new heading.
If you are adding a new heading, please add it at the bottom, but above ``.. _custom:``.

.. code-block:: rst

    Radiation Forces
    ^^^^^^^^^^^^^^^^

    .. _radiation_forces:

    radiation_forces
    ****************

The line of ``^`` characters creates a new subsection, so we're making a subsection named `Radiation Forces` to hold all available radiation force implementations.

Underneath we add the name of our new effect, with a line of ``*`` underneath to create a sub-sub-section.
This should be written lowercase such that we can substitute that effect_name into ``rebx_add_effect_name`` etc. and call the right function! (here rebx_add_radiation_forces).

The ``.. _radiation_forces:`` creates a target that can be used to crossreference to your implementation from other parts of the documentation.

You should now copy paste another documentation entry (e.g., radiation_forces) to make sure you keep the same format.
In our case

.. code-block:: rst

    ======================= ===============================================
    Authors                 H. Rein, D. Tamayo
    Implementation Paper    *In progress*
    Based on                `Burns et al. 1979 <http://labs.adsabs.harvard.edu/adsabs/abs/1979Icar...40....1B/>`_.
    C Example               :ref:`c_example_rad_forces_debris_disk`, :ref:`c_example_rad_forces_circumplanetary`.
    Python Example          `Radiation_Forces_Debris_Disk.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Radiation_Forces_Debris_Disk.ipynb>`_,
                            `Radiation_Forces_Circumplanetary_Dust.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Radiation_Forces_Circumplanetary_Dust.ipynb>`_.
    ======================= ===============================================

This creates a pretty table in the online documentation.
``Authors`` says who wrote the code.
``Implementation paper`` is the paper that you'd like to be cited by people using your implementation.
``Based on`` is the paper that the equations you used come from.

``C Example`` is a link to the C Example you wrote.
All C examples in the ``reboundx/examples`` directory are automatically built into the documentation, and have cross-reference targets of the form ``c_example_foldername``, where foldername is the name of your example folder in ``reboundx/examples``. 

For the ``Python Example`` line, edit the link from another documentation entry with the name of your ipython notebook filename (in both the title and bracketed URL).

Underneath your table, provide a description that will inform users when it's appropriate to apply your effect (and when it's not!).

Finally, if your effect uses a parameters structure and/or particle parameters, we add tables:

.. code-block:: rst

    This applies radiation forces to particles in the simulation.  
    It incorporates both radiation pressure and Poynting-Robertson drag.
    Only particles whose `beta` parameter is set will feel the radiation.  

    **Effect Structure**: *rebx_params_radiation_forces*

    =========================== ==================================================================
    Field (C type)              Description
    =========================== ==================================================================
    c (double)                  Speed of light in the units used for the simulation.
    source_index (int)          Index in the `particles` array for the radiation source.
    =========================== ==================================================================

    **Particle Parameters**

    Only particles with their ``beta`` parameter set will feel radiation forces.

    =========================== ======================================================
    Name (C type)               Description
    =========================== ======================================================
    beta (double)               Ratio of the radiation force to the gravitational force
                                from the radiation source.
    =========================== ======================================================

These are provided as a quick reference for the user.
Replace the tables with ``*None*`` if your effect has no effect structure or no associated particle parameters.

You can check how everything looks by navigating to ``reboundx/doc`` and typing ``make clean``, then ``make html``.
Then navigate to ``reboundx/doc/_build/html`` and open ``index.html`` in your browser.
The main effects page (with the tables) is on the left: REBx Effects & Parameters.
The automatically included documentation will be under API Documentation (Python) and API Documentation (C).

You're Done!
------------

Send me a pull request (:ref:`pullrequest`), or e-mail me your code, and I'd be happy to incorporate it!

