.. _add_effect:

Adding A New Effect
============================

REBOUNDx aims to be modular, so that when adding your own new effect, you don't have to get into all the low-level code details.
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

Writing the C Code
------------------ 

Good documentation is important for new users to actually use your effect. 
REBOUNDx tries to make this process easier by having you document as you write the code, and automating its inclusion into the online documentation.

reboundx.h
^^^^^^^^^^

The place to start is ``reboundx/src/reboundx.h``, which defines the REBOUNDx API (all the functions/data structures the user will use).

First, under ``Parameter structures for each effect``, you should define your effect structure, following the naming convention ``rebx_params_effect``.  In our case

.. code-block:: c

    struct rebx_params_radiation_forces {
         int source_index;                   ///< Index of particle in particles array that is the source of the radiation.
         double c;                           ///< Speed of light in units appropriate for sim->G and initial conditions.
    }

Make sure to use the comment code as above, and it will automatically get built into the documentation!
Then under ``Functions for adding effects``, add the prototype for your effect adder function, which should follow the naming convention ``rebx_add_effect``, should take any effect parameters you want to take from the user, and return a pointer to the structure you defined above.  
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

In ``effect.h``, you only have to modify the file, brief, and author fields at the top, the include guards (ifndef, define lines) and substitute the name of your effect in the function name.  
Everything else should be kept exactly the same.

In ``effect.c``, we first write the effect adder function.
In our case

.. code-block:: c

    struct rebx_params_radiation_forces* rebx_add_radiation_forces(struct rebx_extras* rebx, int source_index, double c){
        struct rebx_params_radiation_forces* params = malloc(sizeof(*params));
        params->c = c;
        params->source_index = source_index;

        sim->force_is_velocity_dependent = 1;
        rebx_add_force(rebx, params, "radiation_forces", rebx_radiation_forces);
        return params;
    }

Is your effect an additional force, or a post timestep modification (i.e., omething to do between REBOUND timesteps)?

In our case (``radiation_forces``), we have an additional force, but for example mass loss would be a post timestep modification. 
