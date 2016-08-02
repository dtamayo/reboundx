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
#ifndef REBX_REBOUNDX_H
#define REBX_REBOUNDX_H

#ifndef M_PI 
#define M_PI 3.1415926535879323846
#endif
#define REBX_C 10064.9150404 // speed of light in default units of AU/(yr/2pi)

#include <stdint.h>
#include "rebound.h"
#include "core.h"
#include "rebxtools.h"

extern const char* rebx_build_str;      ///< Date and time build string.
extern const char* rebx_version_str;    ///<Version string.

/******************************************
  Enums that might be shared across effects
*******************************************/

/**
 * @brief Enumeration for different coordinate systems.
 */
enum REBX_COORDINATES{
    REBX_COORDINATES_JACOBI,                        ///< Jacobi coordinates.  Default for REBOUND/REBOUNDx.
    REBX_COORDINATES_BARYCENTRIC,                   ///< Coordinates referenced to the center of mass of the whole system.
    REBX_COORDINATES_PARTICLE,                      ///< Coordinates referenced to a particular particle.
};

/****************************************
  General REBOUNDx Functions
*****************************************/
/**
 * \name Main REBOUNDx Functions
 * @{
 */
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
 * @brief Disattaches REBOUNDx instance from simulation and resets simulation's function pointers.
 * @param sim Pointer to the simulation from which to remove REBOUNDx
 */
void rebx_remove_from_simulation(struct reb_simulation* sim);

/**
 * @brief Frees all memory allocated by REBOUNDx instance.
 * @details Should be called after simulation is done if memory is a concern.
 * @param rebx the rebx_extras pointer returned from the initial call to rebx_init.
 */
void rebx_free(struct rebx_extras* rebx);
/** @} */
/** @} */

/****************************************
  Functions for adding effects
*****************************************/
/**
 * \name Functions for adding effects in REBOUNDx
 * @{
 */
/**
 * @defgroup EffectAdders
 * @details These are the functions for adding effects in REBOUNDx.
 * @{
 */

/**
 * @brief Main function for adding effects in REBOUNDx.
 * @param rebx Pointer to the rebx_extras instance returned by rebx_init.
 * @param name Name of the effect we want to add.
 * @return Returns a pointer to a rebx_effect structure for the effect.
 */
struct rebx_effect* rebx_add(struct rebx_extras* rebx, const char* name);

/**
 * @brief Function for adding a custom force in REBOUNDx.
 * @param rebx Pointer to the rebx_extras instance returned by rebx_init.
 * @param name String with the name of the custom effect.
 * @param custom_force User-implemented function that updates the accelerations of particles.
 * @param force_is_velocity_dependent Should be set to 1 if the custom force uses particle velocities, 0 otherwise.
 * @return Returns a pointer to a rebx_effect structure for the effect.
 */
struct rebx_effect* rebx_add_custom_force(struct rebx_extras* rebx, const char* name, void (*custom_force)(struct reb_simulation* const sim, struct rebx_effect* const effect), const int force_is_velocity_dependent);

/**
 * @brief Function for adding a custom post_timestep_modification in REBOUNDx.
 * @param rebx Pointer to the rebx_extras instance returned by rebx_init.
 * @param name String with the name of the custom effect.
 * @param custom_ptm User-implemented function that updates particles.
 * @return Returns a pointer to a rebx_effect structure for the effect.
 */
struct rebx_effect* rebx_add_custom_post_timestep_modification(struct rebx_extras* rebx, const char* name, void (*custom_ptm)(struct reb_simulation* const sim, struct rebx_effect* const effect));
/** @} */
/** @} */

/********************************
 * Parameter manipulation functions
 *******************************/

/**
 * \name Functions for accessing and modifying particle and effect parameters.
 * @{
 */
/**
 * @defgroup ParameterManipulators
 * @brief Functions for accessing and modifying particle and effect parameters.
 * @{
 */

/**
 * @brief Removes a parameter from either a reb_particle or rebx_effect structure.
 * @param object Pointer to either a reb_particle or rebx_effect to which to add the parameter.
 * @param param_name Name of the parameter we want to remove.
 * @return 1 if parameter found and successfully removed, 0 otherwise.
 */
int rebx_remove_param(const void* const object, const char* const param_name);

/**
 * @brief Sets a parameter of type double for a parameter in a linked list of rebx_params.
 * @param object Pointer to either a reb_particle or rebx_effect to which to add the parameter.
 * @param param_name Name of the parameter we want to set (see Effects page at http://reboundx.readthedocs.org for what parameters are needed for each effect)
 * @param value Value to which we want to set the parameter.
 */
void rebx_set_param_double(void* object, const char* const param_name, double value);

/**
 * @brief Gets a parameter value of type double from a rebx_param linked list.
 * @param object Pointer to either a reb_particle or rebx_effect to which to add the parameter.
 * @param param_name Name of the parameter we want to get (see Effects page at http://reboundx.readthedocs.org)
 * @return Pointer to the parameter. NULL if parameter is not found in object (user must check for NULL to avoid segmentation fault).
 */
double* rebx_get_param_double(const void* const object, const char* const param_name);

int* rebx_get_param_int(const void* const object, const char* const param_name);
void rebx_set_param_int(void* object, const char* const param_name, int value);

/** @} */
/** @} */

/******************************************
  Convenience functions for various effects
*******************************************/

/**
 * \name Convenience Functions for Various Effects
 * @{
 */
/**
 * @defgroup ConvFunc
 * @brief Convenience functions for different REBOUNDx effects.
 * @{
 */

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
/**
 * @brief Calculates the particle radius from physical parameters and beta, the ratio of radiation to gravitational forces from the star.
 * @param G Gravitational constant.
 * @param c Speed of light.
 * @param source_mass Mass of the source body.
 * @param source_luminosity Luminosity of radiation source.
 * @param beta ratio of radiation force to gravitational force from the radiation source body.
 * @param density density of particle.
 * @param Q_pr Radiation pressure coefficient (Burns et al. 1979).
 * @return Particle radius (double).
 */
double rebx_rad_calc_particle_radius(const double G, const double c, const double source_mass, const double source_luminosity, const double beta, const double density, const double Q_pr);

/**
 * @brief Calculates the hamiltonian for gr_potential, including the classical Hamiltonian.
 * @param sim pointer to the REBOUND simulation
 * @param gr_potential Effect structure returned by rebx_add("gr_potential").
 * @return Total Hamiltonian, including classical Hamiltonian (double).
 */
double rebx_gr_potential_hamiltonian(struct reb_simulation* const sim, const struct rebx_effect* const gr_potential);

/**
 * @brief Calculates the hamiltonian for gr, including the classical Hamiltonian.
 * @details Assumes there is only one source particle (with gr_source set to 1)
 * @param sim pointer to the REBOUND simulation
 * @param gr Effect structure returned by rebx_add("gr").
 * @return Total Hamiltonian, including classical Hamiltonian (double).
 */
double rebx_gr_hamiltonian(struct reb_simulation* const sim, const struct rebx_effect* const gr);

/**
 * @brief Calculates the hamiltonian for gr_full, including the classical Hamiltonian.
 * @param sim pointer to the REBOUND simulation
 * @param gr_full Effect structure returned by rebx_add("gr_full")
 * @return Total Hamiltonian, including classical Hamiltonian (double).
 */
double rebx_gr_full_hamiltonian(struct reb_simulation* const sim, const struct rebx_effect* const gr_full);

/** @} */
/** @} */

#endif
