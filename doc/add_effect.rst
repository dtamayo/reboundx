.. _add_effect:

Adding A New Effect
============================

Whether you're playing around with new physics in your integrations or want to contribute a new effect to REBOUNDx, the process is easy.
Our hope is that as people use the package and work on new problems, they will contribute their new effects so others can use (and cite!) their implementations.
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

Note: Different implementations for the same effect should be added as separate effects.
For example, there are separate files for ``gr``, ``gr_potential``, and ``gr_full``.

Writing the C Code
------------------ 

effect.c
^^^^^^^^^^^^^^^^^^^^^

First copy an existing ``.c`` file in ``reboundx/src``, so that you have the license and right code structure to work from.
Typically you'll be able to reuse many parts of the code.

The function prototype should always stay the same, so we just change the name

.. code-block:: c

    void rebx_stark_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){

All forces are always passed a pointer to the simulation ``sim``, a force structure ``force``, a ``particles`` array, and the number ``N`` of particles in the array (real particles, not counting variational particles--you don't have to worry if you don't know what those are).
The only thing our function needs to do is evaluate the accelerations for each particle, and add those to each particles' acceleration vector, as  we'll do below.
Note that we should specifically update the passed ``particles`` array (NOT ``sim->particles``).

As a simple first example, let's hardcode a constant acceleration along the x direction for the first particle, following the ipynb example:

.. code-block:: c

    void rebx_stark_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
        particles[1].ax += 0.01;
    }

Note that we ADD (+=) to the particles' accelerations, rather than overwrite them (i.e., particles[1].ax = 0.01). This way the various accelerations acting on particles can be accumulated.

core.c and core.h
^^^^^^^^^^^^^^^^^

You need to add your new force as a new ``else if`` in the ``rebx_load_force`` function in reboundx/src/core.c, referencing the function you've written, and the type of force.
If evaluation of your accelerations involves the particle velocities, set ``REBX_FORCE_VEL``, otherwise ``REBX_FORCE_POS``:

.. code-block:: c
    
    ...
    else if (strcmp(name, "stark_force") == 0){
        force->update_accelerations = rebx_stark_force;;
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
^^^^^^^^^

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

        struct reb_particle planet = {0};  # add a planet on a circular orbit (with default units where G=1)
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
(See Sec. 2.4 of `OpenGL Keyboard Commands <http://rebound.readthedocs.org/en/latest/c_quickstart.html>`_ for a list of the visualization keyboard commands).
If you get an error about OpenGL or GLUT, just google `install openGL glut libraries <your OS here>` for instructions, or open your ``Makefile`` and set OPENGL=0 to turn it off.

Python Example
^^^^^^^^^^^^^^

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
If you want to make your effect more flexible, so that users can change parameters at runtime, check out :ref:`advanced_add_effect`, and :ref:`contribute` if you want to add your effect to REBOUNDx so others can also use it.

Operators
^^^^^^^^^

While the above example shows how to add a new force, adding operators is very analogous.
As opposed to updating accelerations, operators should update the particle states (typically their velocities).
Operators make up splitting schemes, and you should read our REBOUNDx paper if you're not familiar with them.

The only difference in implementation from the above is that you would update ``rebx_load_operator`` instead of ``rebx_load_force``, and your function prototype should look like

.. code-block:: c

    void rebx_my_operator(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){

where ``sim`` is again a pointer to the simulation, ``operator`` is an operator struct analogous to the ``force`` struct, and ``dt`` is the length of time over which the operator should act. See ``modify_mass.c`` for an example.
