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

void rebx_modify_orbits_forces(struct reb_simulation* const sim);
void rebx_modify_orbits_direct(struct reb_simulation* const sim);
void rebx_gr_full(struct reb_simulation* const sim);
void rebx_gr_potential(struct reb_simulation* const sim);
void rebx_gr(struct reb_simulation* const sim);
void rebx_radiation_forces(struct reb_simulation* const sim);

extern const char* rebx_build_str;      ///< Date and time build string.
extern const char* rebx_version_str;    ///<Version string.

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
 * @brief Adds post-Newtonian corrections arising from all bodies in the simulation.  Safe, but slow.
 */
void rebx_add_gr(struct rebx_extras* rebx, double c);
/**
 * @brief Adds post-Newtonian corrections arising only from particles[0].  Gets precessions and mean motions right.  Accurate and fast.
 */
void rebx_add_gr_single_mass(struct rebx_extras* rebx, double c);
/**
 * @brief Adds simple potential for post-Newtonian corrections arising only from particles[0].  Gets precessions but not mean motions correct.  Fastest.
 */
void rebx_add_gr_potential(struct rebx_extras* rebx, double c);
/**
 * @brief Adds radiation forces to the simulation (i.e., radiation pressure and Poynting-Robertson drag).
 * @param source Pointer to the particle that is the source of the radiation.
 * @param c Speed of light.
 */
void rebx_add_radiation_forces(struct rebx_extras* rebx, struct reb_particle* source, double c);
/**
 * @brief Allows user to specify their own post timestep modifications. Behaviour is identical to setting up post timestep modifications in REBOUND itself.
 * @param custom_ptm Custom post-timestep modification function.
 */
void rebx_add_custom_ptm(struct rebx_extras* rebx, void (*custom_ptm)(struct reb_simulation* const r) );

/**
 * @brief Allows user to specify their own extra forces. Behaviour is identical to setting up extra forces  in REBOUND itself.
 * @param force_is_velocity_dependent Set to 1 if force is velocity dependent.
 */
void rebx_add_custom_forces(struct rebx_extras* rebx, void (*custom_forces)(struct reb_simulation* const r), int force_is_velocity_dependent);
/** @} */

/**
 * @defgroup GettersSetters
 * @{
 */

// Getter setter landmark for add_param.py
void rebx_set_tau_a(struct reb_particle* p, double value); 
double rebx_get_tau_a(struct reb_particle* p);
void rebx_set_tau_e(struct reb_particle* p, double value); 
void rebx_set_tau_inc(struct reb_particle* p, double value); 
void rebx_set_tau_omega(struct reb_particle* p, double value); 
void rebx_set_tau_Omega(struct reb_particle* p, double value); 
double rebx_get_tau_e(struct reb_particle* p);
double rebx_get_tau_inc(struct reb_particle* p);
double rebx_get_tau_omega(struct reb_particle* p);
double rebx_get_tau_Omega(struct reb_particle* p);
void rebx_set_beta(struct reb_particle* p, double value);
double rebx_get_beta(struct reb_particle* p);

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
 * Internal Functions and data structures
 ***************************************/

/****************************************
Basic types in REBOUNDx
*****************************************/

/*  Enumeration for the different groups of parameters that can be added to particles.*/
enum REBX_PARAMS{
    ORB_TAU,                            // Parameter holding timescales for orbit modifications (migration etc.).
    RAD_BETA,                           // Ratio of radiation pressure force to gravitational force from the star.
};

/*  Main structure used for all parameters added to particles.
    These get added as nodes to a linked list for each particle, stored at particles[i].ap.*/
struct rebx_param{
    void* paramPtr;                     // Pointer to the parameter (void* so it can point to different types of structs).
    enum REBX_PARAMS param_type;        // Identifier for the type of parameter.
    struct rebx_param* next;            // Pointer to the next parameter in the linked list.
};

/*  Nodes for a linked list to all the parameters that have been allocated by REBOUNDx (so it can later free them).*/
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
    double c;                           // Speed of light in units appropriate for sim->G and initial conditions.
};

/*  Structure for adding radiation forces to the simulation.*/
struct rebx_params_radiation_forces{
    struct reb_particle* source;        // Pointer to the particle that's the source of the radiation.
    double c;                           // Speed of light
};

/*************************************************
Structures for effect-specific particle parameters
**************************************************/
/*  Structure to hold all the orbit modification timescales (used by both modify_orbits_direct and modify_orbits_forces).*/
struct rebx_orb_tau{
    double tau_a;                       // Semimajor axis e-folding timescale (<0 = damp, >0 = grow).
    double tau_e;                       // Eccentricity e-folding timescale (<0 = damp, >0 = grow).
    double tau_inc;                     // Inclination e-folding timescale (<0 = damp, >0 = grow).
    double tau_omega;                   // Pericenter precession timescale (linear) (>0 = prograde, <0 = retrograde).
    double tau_Omega;                   // Nodal precession timescale (linear) (>0 = prograde, <0 = retrograde).
};

/****************************************
Main REBOUNDx structure
*****************************************/
struct rebx_extras {    
    struct reb_simulation* sim;                             // Pointer to the simulation REBOUNDx is linked to.
    struct rebx_param_to_be_freed* params_to_be_freed;      // Linked list with pointers to all parameters allocated by REBOUNDx (for later freeing).

    void (**ptm) (struct reb_simulation* const sim);        // Pointer to an array of function pointers for all the post-timestep modifications added by user.
    void (**forces) (struct reb_simulation* const sim);     // Pointer to an array of function pointers for all the force-based modifications added by user.

    int Nptm;                                               // Number of post-timestep modifications currently in simulation.
    int Nforces;                                            // Number of force-based modifications currently in simulation.

    struct rebx_params_modify_orbits modify_orbits_forces;  // Structure for migration/ecc,inc damping/precession using forces.
    struct rebx_params_modify_orbits modify_orbits_direct;  // Structure for migration/ecc,inc damping/precession directly altering orbital elements.
    struct rebx_params_gr gr;                               // Structure for post-Newtonian corrections.
    struct rebx_params_radiation_forces radiation_forces;   // Structure for radiation forces.
};

/****************************************
Internal functions
*****************************************/

/* Main routines called each timestep. */
void rebx_forces(struct reb_simulation* sim);               // Calls all the forces that have been added to the simulation.
void rebx_ptm(struct reb_simulation* sim);                  // Calls all the post-timestep modifications that have been added to the simulation.

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx); // Initializes all pointers and values.

/* Garbage collection routines. */
void rebx_free_params(struct rebx_extras* rebx);            // Steps through linked list to free all allocated parameters.
void rebx_free_pointers(struct rebx_extras* rebx);          // Frees all the remaining pointers in rebx_extras.

/* Internal utility functions. */
void* rebx_search_param(const struct reb_particle* p, enum REBX_PARAMS param);  // returns rebx_param corresponding to the passed param in the passed particle.  If it doesn't exist, returns NULL.
void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param); // add a node for param in the rebx_params_to_be_freed linked list.

/* Internal parameter adders (need a different one for each REBX_PARAM type). */
void rebx_add_param_double(struct reb_particle* p, enum REBX_PARAMS param_type, double value);
void rebx_add_param_orb_tau(struct reb_particle* p);        // add a rebx_orb_tau parameter to particle p

/** @} */

/* Function for testing whether REBOUNDx can load librebound.so and call REBOUND functions. */
double install_test(void);

#endif

