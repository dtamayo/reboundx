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

#ifndef LIBREBX_H
#define LIBREBX_H
#ifndef M_PI 
#define M_PI 3.1415926535879323846
#endif
#define C_DEFAULT 10064.9150404 // speed of light in default units of AU/(yr/2pi)

#include <stdint.h>
#include "rebound.h"
#include "core.h"
#include "rebxtools.h"
#include "modify_orbits_direct.h"

extern const char* rebx_build_str;      ///< Date and time build string.
extern const char* rebx_version_str;    ///<Version string.

/*  Enumeration for different types of coordinate systems.*/
enum REBX_COORDINATES{
    JACOBI,                             // Default.  Uses Jacobi coordinates.
    BARYCENTRIC,                        // Uses coordinates referenced to the center of mass of the whole system.
    HELIOCENTRIC                        // Uses coordinates referenced to sim->particles[0].
};

struct rebx_params_modify_orbits_direct{
    double p;                           // p parameter from Deck & Batygin (2015) for how e-damping couples to a-damping at order e^2.  p=0 : no damping (default), p=1 : e-damping at constant angular momentum.
    enum REBX_COORDINATES coordinates;  // Identifier for the coordinate system that should be used for the damping.
};

struct rebx_params_modify_orbits_forces{
    enum REBX_COORDINATES coordinates;  // Identifier for the coordinate system that should be used for the damping.
};

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
 * @brief Adds orbit modifications (migration, damping, precession), altering the orbital elements directly.
 */
struct rebx_params_modify_orbits_direct* rebx_add_modify_orbits_direct(struct rebx_extras* rebx);
/**
 * @brief Adds orbit modifications (migration, damping, precession), implemented as forces.
 */
struct rebx_params_modify_orbits_forces* rebx_add_modify_orbits_forces(struct rebx_extras* rebx);
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
void rebx_gr_full(struct reb_simulation* const sim, struct rebx_effect* const gr);
void rebx_gr_potential(struct reb_simulation* const sim, struct rebx_effect* const gr);
void rebx_gr(struct reb_simulation* const sim, struct rebx_effect* const gr);
void rebx_radiation_forces(struct reb_simulation* const sim, struct rebx_effect* const rad);


#endif
