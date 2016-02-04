/**
 * @file    reboundx.h
 * @brief   REBOUNDx API definition.
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _LIBREBX_H
#define _LIBREBX_H
#ifndef M_PI 
#define M_PI 3.1415926535879323846
#endif
#define C_DEFAULT 10064.9150404 // speed of light in default units of AU/(yr/2pi)

#include "rebound.h"
#include "rebxtools.h"


extern const char* rebx_build_str;      ///< Date and time build string.
extern const char* rebx_version_str;    ///<Version string.
/****************************************
 * Internal Functions and data structures
 ***************************************/

/****************************************
Basic types in REBOUNDx
*****************************************/

/*  Enumeration for the different groups of parameters that can be added to particles.*/
enum REBX_PARAMS{
    RAD_BETA,                               // Ratio of radiation to gravitational force (Burns et al. 1979)
    TAU_LITTLE_OMEGA,                   // Period of linear apsidal precession/regression
    TAU_BIG_OMEGA,                      // Period of linear nodal precession/regression
    TAU_INC,                            // Inclination exponential growth/damping timescale
    TAU_E,                              // Eccentricity exponential growth/damping timescale
    TAU_A,                              // Semimajor axis exponential growth/damping timescale
};

/*	Enumeration for the different effects that can be added in REBOUNDx.  See reboundx.readthedocs.org for details on the implementation.*/
enum REBX_EFFECTS{
	RADIATION_FORCES,
	MODIFY_ORBITS_DIRECT,
	MODIFY_ORBITS_FORCES,
	GR,
	GR_FULL,
	GR_POTENTIAL,
    CUSTOM_POST_TIMESTEP_MODIFICATION,
};

/* 	Main structure used for all parameters added to particles.
 	These get added as nodes to a linked list for each particle, stored at particles[i].ap.*/
struct rebx_param{
    void* paramPtr;                     // Pointer to the parameter (void* so it can point to different types of structs).
    enum REBX_PARAMS param_type;        // Identifier for the type of parameter.
    struct rebx_param* next;            // Pointer to the next parameter in the linked list.
};

/*  Structure for all REBOUNDx effects.
 *  These get added as nodes to the effects linked list in the rebx_extras structure.*/
struct rebx_effect{
	void* paramsPtr;						// Pointer to the effect params structure (void* so it can point to different effect structs).
	void (*functionPtr) (struct reb_simulation* const sim, struct rebx_effect* const effect);	// Pointer to the function to carry out the additional effect.
	enum REBX_EFFECTS effect_type;		// Identifier for the type of effect.
	struct rebx_effect* next;			// Pointer to the next effect in the linked list.
};

/*	Nodes for a linked list to all the parameters that have been allocated by REBOUNDx (so it can later free them).*/
struct rebx_param_to_be_freed{
    struct rebx_param* param;           // Pointer to a parameter node allocated by REBOUNDx.
    struct rebx_param_to_be_freed* next;// Pointer to the next node in the linked list rebx_extras.params_to_be_freed.
};

/****************************************
Enums and structs for the particular modifications.
*****************************************/

/*  Enumeration for different types of coordinate systems.*/
enum REBX_COORDINATES{
    JACOBI,                             // Default.  Uses Jacobi coordinates.
    BARYCENTRIC,                        // Uses coordinates referenced to the center of mass of the whole system.
    HELIOCENTRIC                        // Uses coordinates referenced to sim->particles[0].
};

/* Structure for orbit modifications (modify_orbits_direct and modify_orbits_forces).*/
struct rebx_params_modify_orbits{
    double p;                           // p parameter from Deck & Batygin (2015) for how e-damping couples to a-damping at order e^2.  p=0 : no damping (default), p=1 : e-damping at constant angular momentum.
    enum REBX_COORDINATES coordinates;  // Identifier for the coordinate system that should be used for the damping.
};

/*  Structure for adding post-Newtonian corrections.*/
struct rebx_params_gr {
    int source_index;                   // Index of particle in particles array causing GR corrections.
    double c;                           // Speed of light in units appropriate for sim->G and initial conditions.
};

/*  Structure for adding radiation forces to the simulation.*/
struct rebx_params_radiation_forces{
    int source_index;                   // Index of particle in particles array that provides radiation.
    double c;                           // Speed of light
};

/****************************************
Main REBOUNDx structure
*****************************************/
struct rebx_extras {	
	struct reb_simulation* sim;								// Pointer to the simulation REBOUNDx is linked to.
	struct rebx_effect* post_timestep_modifications;		// Linked list with pointers to all the post-timestep modifications added to the simulation.
	struct rebx_effect* forces;                             // Linked list with pointers to all the additional forces added to the simulation.
	struct rebx_param_to_be_freed* params_to_be_freed; 		// Linked list with pointers to all parameters allocated by REBOUNDx (for later freeing).

};

/****************************************
Internal functions
*****************************************/

/* Main routines called each timestep. */
void rebx_forces(struct reb_simulation* sim);               // Calls all the forces that have been added to the simulation.
void rebx_post_timestep_modifications(struct reb_simulation* sim);                  // Calls all the post-timestep modifications that have been added to the simulation.

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx); // Initializes all pointers and values.

/* Garbage collection routines. */
void rebx_free_params(struct rebx_extras* rebx);            // Steps through linked list to free all allocated parameters.
void rebx_free_effects(struct rebx_effect* effects);        // Frees all effects in effects linked list (forces or post_timestep_modifications)
void rebx_free_pointers(struct rebx_extras* rebx);          // Frees all the remaining pointers in rebx_extras.

/* Internal utility functions. */
void* rebx_get_param(const struct reb_particle* p, enum REBX_PARAMS param);  // returns rebx_param corresponding to the passed param in the passed particle.  If it doesn't exist, returns NULL.
struct rebx_effect* rebx_get_effect_in(struct rebx_effect* effects, enum REBX_EFFECTS effect_type);
// returns the first effect in the linked list effects that matches effect_type
struct rebx_effect* rebx_get_effect(struct rebx_extras* rebx, enum REBX_EFFECTS effect_type);   // searches both forces and post_timestep_modifications and returns first effect mathing effect_type
void* rebx_get_effect_params_in(struct rebx_effect* effects, enum REBX_EFFECTS effect_type);
// returns the pointer to the parameters for the first effect in the linked list effects that matches effect_type
    
void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param); // add a node for param in the rebx_params_to_be_freed linked list.

/* Internal parameter adders (need a different one for each REBX_PARAM type). */
void rebx_add_param_double(struct reb_particle* p, enum REBX_PARAMS param_type, double value);

/** @} */

/* Function for testing whether REBOUNDx can load librebound.so and call REBOUND functions. */
double install_test(void);

/****************************************
  User API
*****************************************/

/**
 * @defgroup MainRebxFunctions 
 * @details These are the top level routines that one needs when using REBOUNDx.
 * @{
 */
/**
 * @brief Adds REBOUNDx functionality to a passed REBOUND simulation.
 * @details Allocates memory for a REBOUNDx structure, initializes all variables and returns a pointer to
 * the rebx_extras structure.  The function must be called before calling any other REBOUNDx functions.
 * @param sim reb_simulation pointer to the simulation that you want to add REBOUNDx functionality.
 * @return Returns a pointer to a rebx_extras structure, which holds all the information REBOUNDx needs.
 */
struct rebx_extras* rebx_init(struct reb_simulation* sim);
/**
 * @brief Frees all memory allocated by REBOUNDx instance.
 * @details Should be called after simulation is done if memory is a concern.
 * @param rebx the rebx_extras pointer returned from the initial call to rebx_init.
 */
void rebx_free(struct rebx_extras* rebx);
/** @} */

/**
 * @defgroup AddEffect
 * @{
 */
/**
 * @brief Adds orbit modifications (migration, damping, precession), implemented as forces.
 */
void rebx_add_modify_orbits_forces(struct rebx_extras* rebx);
/**
 * @brief Adds orbit modifications (migration, damping, precession), altering the orbital elements directly.
 */
void rebx_add_modify_orbits_direct(struct rebx_extras* rebx);
/**
 * @brief Adds post-Newtonian corrections arising only from a single particle (specified by source).  Gets precessions and mean motions right.  Accurate and fast.
 */
void rebx_add_gr(struct rebx_extras* rebx, struct reb_particle* source, double c);
/**
 * @brief Adds post-Newtonian corrections arising from all bodies in the simulation.  Safe, but slower.
 */
void rebx_add_gr_full(struct rebx_extras* rebx, double c);
/**
 * @brief Adds simple potential for post-Newtonian corrections arising only from a single particle (specified by source).  Gets precessions but not mean motions correct.  Fastest.
 */
void rebx_add_gr_potential(struct rebx_extras* rebx, struct reb_particle* source, double c);
/**
 * @brief Adds radiation forces to the simulation (i.e., radiation pressure and Poynting-Robertson drag).
 * @param source particle that is the source of the radiation.
 * @param c Speed of light.
 */
void rebx_add_radiation_forces(struct rebx_extras* rebx, struct reb_particle* source, double c);
/**
 * @brief Allows user to specify their own post timestep modifications. Behaviour is identical to setting up post timestep modifications in REBOUND itself.
 * @param custom_post_timestep_modifications Custom post-timestep modification function.
 * @param custom_params Custom parameters container.  Pass NULL if you don't want to use it.
 */
void rebx_add_custom_post_timestep_modification(struct rebx_extras* rebx, void (*custom_post_timestep_modification)(struct reb_simulation* const sim, struct rebx_effect* const custom_effect), void* custom_params);

/**
 * @brief Allows user to specify their own extra forces. Behaviour is identical to setting up extra forces  in REBOUND itself.
 * @param custom_forces Custom forces function.
 * @param force_is_velocity_dependent Set to 1 if force is velocity dependent.
 * @param custom_params Custom parameters container.  Pass NULL if you don't want to use it.
 */
void rebx_add_custom_force(struct rebx_extras* rebx, void (*custom_force)(struct reb_simulation* const sim, struct rebx_effect* const custom_effect), int force_is_velocity_dependent, void* custom_params);
/** @} */

/**
 * @defgroup GettersSette
 * @{
 */

// Getter setter landmark for add_param.py
void rebx_set_beta(struct reb_particle* p, double value);
double rebx_get_beta(struct reb_particle* p);
void rebx_set_tau_omega(struct reb_particle* p, double value);
double rebx_get_tau_omega(struct reb_particle* p);
void rebx_set_tau_Omega(struct reb_particle* p, double value);
double rebx_get_tau_Omega(struct reb_particle* p);
void rebx_set_tau_inc(struct reb_particle* p, double value);
double rebx_get_tau_inc(struct reb_particle* p);
void rebx_set_tau_e(struct reb_particle* p, double value);
double rebx_get_tau_e(struct reb_particle* p);
void rebx_set_tau_a(struct reb_particle* p, double value);
double rebx_get_tau_a(struct reb_particle* p);

/**
 * @brief Searches all added effects and returns the parameters for the FIRST encountered effect of effect_type
 * @param effect_type enum for the type of effect to search for (take name from call to rebx_add, in all caps, e.g., GR_FULL)
 */
void* rebx_get_effect_params(struct rebx_extras* rebx, enum REBX_EFFECTS effect_type);

/** @} */

/**
 * @defgroup ConvFunc
 * @{
 */

/**
 * @brief Calculates beta, the ratio between the radiation pressure force and the gravitational force from the star.
 */
double rebx_rad_calc_beta(struct rebx_extras* rebx, double particle_radius, double density, double Q_pr, double L);
/**
 * @brief Calculates the particle radius from physical parameters and beta, the ratio of radiation to gravitational forces from the star.
 */
double rebx_rad_calc_particle_radius(struct rebx_extras* rebx, double beta, double density, double Q_pr, double L);

/** @} */

/****************************************
  Function prototypes
*****************************************/

void rebx_modify_orbits_forces(struct reb_simulation* const sim, struct rebx_effect* const mod);
void rebx_modify_orbits_direct(struct reb_simulation* const sim, struct rebx_effect* const mod);
void rebx_gr_full(struct reb_simulation* const sim, struct rebx_effect* const gr);
void rebx_gr_potential(struct reb_simulation* const sim, struct rebx_effect* const gr);
void rebx_gr(struct reb_simulation* const sim, struct rebx_effect* const gr);
void rebx_radiation_forces(struct reb_simulation* const sim, struct rebx_effect* const rad);

#endif
