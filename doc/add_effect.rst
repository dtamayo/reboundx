.. _add_effect:

Adding A New Effect
===================

Whether you're playing around with new physics in your integrations or want to contribute a new effect to REBOUNDx, the process is easy.

The :ref:`basic_force` section shows you how to get your new effect working in C and python within 10 minutes. 
The :ref:`adding_parameters` section shows you how to allow the user to set effect and particle parameters at runtime, and the :ref:`contributing` section goes over how to add your new effect to the REBOUNDx repository, so others can use it and find it in the documentation.

Our hope is that as people use the package and work on new problems, they will contribute their new effects, so others can use (and cite!) their implementations.
Contributing the effects you implement (and may want to use in the future) to REBOUNDx has the added benefit of ensuring that they stay up to date as REBOUND and REBOUNDx expand.

*Do I Have To Write It In C?*

Forces and operators are called every timestep, and the overhead of REBOUND calling a Python function each timestep makes it a factor of a few slower than if the effect was written in C.
Therefore all effects in REBOUNDx are written in C.
Often you might want to quickly try something out in Python, and the Custom_Effects.ipynb example shows you how to do that.
Here we walk through an example of porting the stark force in Custom_Effect.ipynb into C.

*Writing Forces*

There are two kinds of effects, forces and operators.
We first consider a force, which boils down to writing a C function that will evaluate and add the relevant accelerations to all the particles.
Behind the scenes, REBOUNDx then takes care of when your acceleration function gets called depending on the REBOUND integrator the user has chosen.
Adding operators is very analogous, and we discuss the differences in :ref:`operators`.

Note: Different implementations for the same effect should be added as separate effects.
For example, there are separate files for ``gr``, ``gr_potential``, and ``gr_full``.

.. _basic_force:

Basic Force
^^^^^^^^^^^

*effect.c*

First copy an existing ``.c`` file in ``reboundx/src``, so that you have the license and right code structure to work from.
Typically, you'll be able to reuse many parts of the code.

The function prototype should always stay the same, so we just change the name

.. code-block:: c

    void rebx_stark_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){

All forces are always passed as a pointer to the simulation ``sim``, a force structure ``force``, a ``particles`` array, and the number ``N`` of particles in the array (real particles, not counting variational particles--you don't have to worry if you don't know what those are).
The only thing our function needs to do is evaluate the accelerations for each particle, and add those to each particles' acceleration vector, as  we'll do below.
Note that we should specifically update the passed ``particles`` array (NOT ``sim->particles``).

As a simple first example, let's hardcode a constant acceleration along the x direction for the first particle, following the ipynb example:

.. code-block:: c

    void rebx_stark_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
        particles[1].ax += 0.01;
    }

Note that we ADD (+=) to the particles' accelerations, rather than overwrite them (i.e., particles[1].ax = 0.01). This way the various accelerations acting on particles can be accumulated.

*core.c and core.h*

You need to add your new force as a new ``else if`` in the ``rebx_load_force`` function in reboundx/src/core.c, referencing the function you've written, and the type of force.
If evaluation of your accelerations involves the particle velocities, set ``REBX_FORCE_VEL``, otherwise ``REBX_FORCE_POS``:

.. code-block:: c
    
    ...
    else if (strcmp(name, "stark_force") == 0){
        force->update_accelerations = rebx_stark_force;
        force->force_type = REBX_FORCE_POS;
    }
    else{
    ...

You also need to add your function prototype at the bottom of reboundx/src/core.h under Force prototypes:

.. code-block:: c

    ...
    void rebx_stark_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N);

Now

.. code-block:: bash

    cd your_path_to/reboundx/scripts
    python add_new_effect.py

This script updates all the makefiles and the pip installation file to include your new effect.

That's it! Your new force now works from both C and Python. Let's try it out.

C Example
*********

Let's test this first in C. This could then turn into a C example for others if you contributed it to REBOUNDx (all REBOUNDx effects have corresponding C examples).
Navigate to the ``reboundx/examples`` folder, and copy any folder to a new one named ``stark_force``.
Now we just modify the ``problem.c`` file in our new ``stark_force`` folder, e.g.:

.. code-block:: c

    #include "rebound.h"
    #include "reboundx.h"
    
    int main(int argc, char* argv[]){
        struct reb_simulation* sim = reb_create_simulation();
        struct reb_particle star = {0}; 
        star.m     = 1.;   
        reb_add(sim, star); 

        struct reb_particle planet = {0};  // add a planet on a circular orbit (with default units where G=1)
        planet.x = 1.;
        planet.vy = 1.;
        reb_add(sim, planet);

        struct rebx_extras* rebx = rebx_attach(sim);  // first initialize rebx
        struct rebx_force* stark = rebx_load_force(rebx, "stark_force"); // add our new force
        rebx_add_force(rebx, stark);                
        
        double tmax = 100000.;
        reb_integrate(sim, tmax);
    }

In the terminal in the ``stark_force`` folder then just ``make clean``, ``make``  and then run it with  ``./rebound``. 
In the visualization press 'w' to see the orbits. You should see a mess with the orbit getting more and less eccentric.
If you get an error about OpenGL or GLUT, just google `install openGL glut libraries <your OS here>` for instructions, or open your ``Makefile`` and set OPENGL=0 to turn it off.

Python Example
**************

Our new effect will now work out of the box without any extra python code.
We just need to make sure that whenever we change C code (like we did above), we reinstall REBOUNDx, i.e. ``pip install -e .`` in the root ``reboundx`` directory. 
Then, e.g. in a jupyter notebook:

..  code-block:: python

    import rebound
    import reboundx
    
    sim = rebound.Simulation()
    sim.add(m=1.)
    sim.add(a=1.)

    rebx = reboundx.Extras(sim)
    stark = rebx.load_force("stark_force")
    rebx.add_force(stark)

    sim.integrate(1e5)

will run with our new effect. We could plot the eccentricity vs time, just like in the Custom_Effects.ipynb ipython_example where we code the effect in python (and is a factor of a few slower than our new C code).

That's all there is to it. 
If you want to make your effect more flexible, so that users can change parameters at runtime, check out :ref:`adding_parameters`, and :ref:`contributing` if you want to add your effect to REBOUNDx so others can also use it.

.. _operators:

Operators
*********

While the above example shows how to add a new force, adding operators is very analogous.
As opposed to updating accelerations, operators should update the particle states (typically their velocities).
Operators make up splitting schemes, and you should read our REBOUNDx paper if you're not familiar with them.

The only difference in implementation from the above is that you would update ``rebx_load_operator`` instead of ``rebx_load_force``, and your function prototype should look like

.. code-block:: c

    void rebx_my_operator(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){

where ``sim`` is again a pointer to the simulation, ``operator`` is an operator struct analogous to the ``force`` struct, and ``dt`` is the length of time over which the operator should act. See ``modify_mass.c`` for an example.

.. _adding_parameters:

Adding Parameters
^^^^^^^^^^^^^^^^^

In :ref:`basic_force` we went over how to add a simple new force, where we simply hardcoded which particles were effected, and by how much. 
REBOUNDx also makes it easy to add parameters to forces and particles so that the user can have the flexibility to choose these values at runtime, can write a script that sets parameters individually on particles, can inspect them to write output to files, change values halfway through, etc.
This will show you how to do that with your effect.

Particularly if you are used to the Python side of REBOUND/REBOUNDx, you should read :ref:`c_quickstart` to see how to access parameters, and it can be very useful to look at forces that are already implemented.

First, you should decide whether force parameters belong on your force or on particles.
For example, the ``radiation_forces`` effect needs to know the speed of light (which will vary if the user changes units), so ``c`` is a parameter that is the same for all particles, and is added to the force.
In our case, if our constant stark acceleration was the same for all particles, we might add it to the force.
If each particle could feel a different acceleration, we would add them to the particles.
That will depend on the physics you're trying to put in--let's add the parameter to the particles as an example.

The first thing to do is register the parameter name in ``rebx_register_default_params`` in ``src/core.c``.
You cannot use particle parameter names that are in use by other effects, so search first for the name you are planning to add.
In order to avoid clashes, we have implemented a convention for new effects that any parameters must start with the acronym for the effect.
So for example the tau parameter for ``tides_constant_time_lag`` is ``tctl_tau``.

You also have to specify the type.
The vast majority of parameters will be integers (REBX_TYPE_INT) or doubles (REBX_TYPE_DOUBLE). 
We go through what to do with new custom types below.

Here let's call our parameter ``stark_acc``, and it should be a double:

.. code-block:: c

    ...
    rebx_register_param(rebx, "stark_acc", REBX_TYPE_DOUBLE);

The user will now be able to set and check the value of this parameter on all particles.
Now we have to do something with it in our ``stark_force`` implementation, following the basic example in :ref:`add_effect`:

.. code-block:: c

    void rebx_stark_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
        struct rebx_extras* const rebx = sim->extras;
        for (int i=0; i<N; i++){
            const double* stark_acc = rebx_get_param(rebx, particles[i].ap, "stark_acc");
            if (stark_acc != NULL){
                particles[i].ax += *stark_acc;
            }
        }
    }

We now iterate through the particle list, check whether each one has its ``stark_acc`` param, and if so, update its x acceleration.
It's good practice to always check the parameter pointers you get back from ``rebx_get_param`` for NULL, since otherwise you will get a seg fault when you dereference them if they have not been set by the user.

C Example
*********

Now in C we can set our particle parameter in our ``problem.c`` file:


.. code-block:: c

    rebx_set_param_double(rebx, &sim->particles[1].ap, "stark_acc", 0.01);

and everything will work as before.

Python Example
**************

Since we've edited the C code, to use it from Python we have to go back to the root ``reboundx`` directory and ``pip install -e .``.
After that, we can access and change our new particle parameters out of the box:

..  code-block:: python

    sim.particles[1].params['stark_acc'] = 0.01

That's it!

Adding Custom Types
*******************

Most of the time, you'll be using integer and double types for the parameters.
But there may be times where you want to, e.g., use a custom struct.
There is a catchall void pointer type (REBX_TYPE_POINTER) for such cases.
This is convenient in C (see the bottom of the ``reboundx/examples/parameters/problem.c`` example), but in python involves casting things manually (see the bottom of ``reboundx/ipython_examples/GettingStartedParameters.ipynb``).

Here we will take that ipython_example with a made up SPH_sim struct and show how to make your new struct easily accessible in python.

First in ``src/reboundx.h``, we need to add a new enum to ``rebx_param_type``:

.. code-block:: c

    enum rebx_param_type{
        REBX_TYPE_NONE,
        ...
        REBX_TYPE_SPHSIM
    };

Then we need to define this struct below under 'Basic types in REBOUNDx':

.. code-block:: c

    struct rebx_SPH_sim {
        double dt;
        int Nparticles;
    };

Then in ``src/core.c``, under ``rebx_register_default_params``, we need to register it with its new type:

.. code-block:: c
    
    void rebx_register_default_params(struct rebx_extras* rebx){
        ...
        rebx_register_param(rebx, "sph_sim", REBX_TYPE_SPHSIM);


On the Python side, at the bottom of ``reboundx/reboundx/extras.py`` we then have to define the ctypes Structure that matches our C structure (google ctypes documentation or follow the existing examples):

.. code-block:: python
    
    class SPH_sim(Structure):
        _fields_ = [("dt",  c_double),
                    ("Nparticles", c_int)]

and on the line below we have to update the mapping ``REBX_C_TO_CTYPES``, which goes from the ``rebx_param_type`` enum (the first thing we edited in this section) to the Python ctypes structure that we just created (SPH_sim). 
The order in this list must match exactly with what's in the ``rebx_param_type`` enum.

.. code-block:: python
    
    REBX_C_TO_CTYPES = [["REBX_TYPE_NONE", None], ["REBX_TYPE_DOUBLE", c_double], ["REBX_TYPE_INT",c_int], ["REBX_TYPE_POINTER", c_void_p], ["REBX_TYPE_FORCE", Force], ["REBX_TYPE_UNIT32", c_uint32], ["REBX_TYPE_ORBIT", rebound.Orbit], ["REBX_TYPE_SPHSIM", SPH_sim]]

Finally, in ``reboundx/reboundx/params.py``, we have to import our new structure and add a matching if clause in ``__setitem__``:

.. code-block:: python

    from .extras import SPH_sim
    ...

    if ctype == SPH_sim:
        if not isinstance(value, SPH_sim):
            raise AttributeError("REBOUNDx Error: Parameter '{0}' must be assigned a SPHsim object.".format(key))
        clibreboundx.rebx_set_param_pointer(self.rebx, byref(self.ap), c_char_p(key.encode('ascii')), byref(value))
    

Now, if we follow the ``ipython_examples/GettingStartedParameters.ipynb``, rather than manually casting things under Custom Parameters, we can simply do

.. code-block:: python

    from reboundx.extras import SPH_sim
    my_sph_sim = SPH_sim()
    my_sph_sim.dt = 0.1
    my_sph_sim.Nparticles = 10000
    gr.params['sph_sim'] = my_sph_sim
    gr.params['sph_sim'].Nparticles

which will output 10000.

Note that these custom structs will still not be written to REBOUNDx binaries.
If this is important to you, feel free to get in touch.

.. _contributing:

Contributing your effect to REBOUNDx
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You've followed the previous steps in :ref:`add_effect` and :ref:`adding_parameters`, and now want to add it to REBOUNDx so others can use (and cite!) your implementation. Great!

A checklist:

* Have you added a C example in ``reboundx/examples/``?
* Have you added a python example in ``reboundx/ipython_examples/``?
* Have you added documentation, so people can find your effect? (see below)
* Is your code machine independent? (see below)

I'm happy to help with any of these. Once you're ready, send me a pull request (see bottom of this page)

Add Your Effect to the Main Documentation Page!
***********************************************

You add the documentation for your effect directly within your ``effect.c`` file.
It will then automatically get built into the :ref:`effects` page.
Easiest is if you copy-paste from another effect source file.

At the top of the comment block, you should edit the file, brief and author lines.
The rest of the documentation goes Below the dollar signs.
For example, for the stark_force implementation we did: 
 
.. code-block:: rst

     * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     *
     * $Dcoumentation Examples$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
     *
     * ======================= ===============================================
     * Authors                 Jane Done
     * Implementation Paper    `Doe and Smith, 2019 <http://labs.adsabs.harvard.edu/adsabs/abs/2019Nature...30...12/>`_,
     * Based on                `Newton and Halley 1692 <http://labs.adsabs.harvard.edu/adsabs/abs/1692/>`_.
     * C Example               :ref:`c_example_stark_force`
     * Python Example          `Stark_Force.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Stark_Force.ipynb>`_,
     * ======================= ===============================================
     * 
     * This applies a constant acceleration along the x direction (Stark force)
     * 
     * **Effect Parameters**
     * 
     * ============================ =========== ==================================================================
     * Field (C type)               Required    Description
     * ============================ =========== ==================================================================
     * None                         -           - 
     * ============================ =========== ==================================================================
     *
     * **Particle Parameters**
     *
     * Any particles with their stark_acc parameter set will feel the corresponding acceleration along x
     *
     * ============================ =========== ==================================================================
     * Field (C type)               Required    Description
     * ============================ =========== ==================================================================
     * stark_acc (double)           No          Size of the acceleration along the x direction
     * ============================ =========== ==================================================================
     * 
     */

We first add the group that our effect belongs to, between dollar signs.
This keeps different implementations of, e.g., general relativity corrections in the same place.
Here I made up a new one called $Documentation Examples$.
If you want to make a new category like here, you have to add it to the /reboundx/doc/effect_headers.rst file.

When you create a new category in that file, you can optionally add a description general to all implementations in the category following the format in the file, which will show up in :ref:`effects`.

.. code-block:: rst

    $$$$$$$$$$$$$$$$$$$$$$
    Documentation Examples 
    ^^^^^^^^^^^^^^^^^^^^^^
    These are effects that have been added as documentation examples

You can also compare with the Orbit Modifications category in that file and how it shows up in the list of effects in the documentation at :ref:`effects`.

Then fill in the table:
``Authors`` says who wrote the code.
``Implementation paper`` is the paper that you'd like to be cited by people using your implementation.
``Based on`` is the paper that the equations you used come from.

``C Example`` is a link to the C Example you wrote.
All C examples in the ``reboundx/examples`` directory are automatically built into the documentation, and have cross-reference targets of the form ``c_example_foldername``, where foldername is the name of your example folder in ``reboundx/examples``. Here it's ``c_example_star_force``.

For the ``Python Example`` line, edit the link from another documentation entry with the name of your ipython notebook filename (in both the title and bracketed URL).

Underneath your table, provide a description that will inform users when it's appropriate to apply your effect (and when it's not!).

Finally, if your effect requires the user to set (possibly optionally) particular effect or particle parameters, we create tables for them too. 

To check how everything looks the way it should, you need to 

.. code-block:: bash

    pip install breathe sphinx

and you need to install `doxygen <http://www.doxygen.nl/manual/install.html>`_. Then

.. code-block:: bash

    cd reboundx/doc/doxygen
    doxygen Doxyfile
    cd reboundx/doc
    make clean
    make html

Then navigate to ``reboundx/doc/_build/html`` and open ``index.html`` in your browser.
The main effects page (with the tables) is on the left: REBx Effects & Parameters.
The automatically included documentation will be under API Documentation (Python) and API Documentation (C).

Is your code machine independent?
*********************************

This is not a requirement, but worth thinking about, given that the rest of REBOUND is machine  independent, allowing anyone to replicate one another's integrations.
In short, the C99 standard guarantees that arithmetic operations (+,-,*,/) and the sqrt function are machine independent. 
All other math library functions (e.g., sin, cos, exp etc.) are heavily optimized for hardware and can give different results (in the last bit) between architectures.
If you can find a way to write your function to only use basic operations, you can be confident that your code is machine independent.

.. _pullrequest:

Putting together a Pull Request
*******************************

If you'd rather e-mail me your code, I'm happy to incorporate it, but if you'd like for GitHub to show your account as a contributor to the project, send me a pull request! 

If you've never done this before, follow the instructions at `Time to Submit Your First PR <http://www.thinkful.com/learn/github-pull-request-tutorial/Expect-a-Thorough-Review#Time-to-Submit-Your-First-PR>`_ up until "Tadaa!" to fork the REBOUNDx repository and make your own local branch.

Now you can modify the code as described below, and can incrementally commit changes.
As a starting point, you can check out `this guide <https://www.atlassian.com/git/tutorials/saving-changes>`_.

After working through this document and making all the changes, you can then send me a pull request by following the rest of the instructions in the pull request tutorial above.
We're always happy to help. Let us know if you have any questions or suggestions for how to improve this tutorial by opening an issue on the REBOUNDx `GitHub page! <https://github.com/dtamayo/reboundx>`_.
