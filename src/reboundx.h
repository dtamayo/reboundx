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
#define C_DEFAULT 10064.9150404 // speed of light in default units of AU/(yr/2pi)

#include <stdint.h>
#include "rebound.h"
#include "core.h"
#include "rebxtools.h"
#include "modify_orbits_direct.h"

extern const char* rebx_build_str;      ///< Date and time build string.
extern const char* rebx_version_str;    ///<Version string.

/******************************************
  Enums that might be shared across effects
*******************************************/

/**
 * @brief Enumeration for different coordinate systems.
 */
enum REBX_COORDINATES{
    JACOBI,                             ///< Jacobi coordinates.  Default for REBOUND/REBOUNDx.
    BARYCENTRIC,                        ///< Coordinates referenced to the center of mass of the whole system.
    HELIOCENTRIC                        ///< Coordinates referenced to particles[0] in the simulation.
};

/*****************************************
  Parameter structures for each effect
******************************************/

struct rebx_params_modify_orbits_direct{
    double p;                           ///< Coupling parameter between eccentricity and semimajor axis evolution.
    enum REBX_COORDINATES coordinates;  ///< Coordinate system that should be used for the calculations.
};

struct rebx_params_modify_orbits_forces{
    enum REBX_COORDINATES coordinates;  ///< Coordinate system that should be used for the calculations.
};

struct rebx_params_gr {
    int source_index;                   ///< Index of particle in particles array causing GR corrections.
    double c;                           ///< Speed of light in units appropriate for sim->G and initial conditions.
};

struct rebx_params_gr_potential {
    int source_index;                   ///< Index of particle in particles array causing GR corrections.
    double c;                           ///< Speed of light in units appropriate for sim->G and initial conditions.
};

struct rebx_params_gr_full {
    double c;                           ///< Speed of light in units appropriate for sim->G and initial conditions.
};

struct rebx_params_radiation_forces {
    int source_index;                   ///< Index of particle in particles array that is the source of the radiation.
    double c;                           ///< Speed of light in units appropriate for sim->G and initial conditions.
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
 * \name Adder Functions for Each REBOUNDx Effect
 * @{
 */
/**
 * @defgroup AddEffect
 * @{
 */
/**
 * @brief Adds orbit modifications, altering the orbital elements directly.
 * @details Silently sets the effect struct's p member to 0, and coordinates to JACOBI.
 * User can modify if needed.
 * @param rebx pointer to the rebx_extras instance
 */
struct rebx_params_modify_orbits_direct* rebx_add_modify_orbits_direct(struct rebx_extras* rebx);

/**
 * @brief Adds orbit modifications, implemented as forces.
 * @details Silently sets the effect struct's coordinates member to JACOBI.
 * User can modify if needed.
 * @param rebx pointer to the rebx_extras instance
 */
struct rebx_params_modify_orbits_forces* rebx_add_modify_orbits_forces(struct rebx_extras* rebx);

/**
 * @brief Adds post-Newtonian corrections arising only from a single particle.  
 * @param rebx pointer to the rebx_extras instance
 * @param source_index Index in the particles array of the body that causes the GR corrections.
 * @param c Speed of light.
 */
struct rebx_params_gr* rebx_add_gr(struct rebx_extras* rebx, int source_index, double c);

/**
 * @brief Adds simple potential for post-Newtonian corrections arising only from a single particle.
 * @param rebx pointer to the rebx_extras instance
 * @param source_index Index in the particles array of the body that causes the GR corrections.
 * @param c Speed of light.
 */
struct rebx_params_gr_potential* rebx_add_gr_potential(struct rebx_extras* rebx, int source_index, double c);

/**
 * @brief Adds post-Newtonian corrections arising from all bodies in the simulation.
 * @param rebx pointer to the rebx_extras instance
 * @param c Speed of light.
 */
struct rebx_params_gr_full* rebx_add_gr_full(struct rebx_extras* rebx, double c);

/**
 * @brief Adds radiation forces to the simulation (i.e., radiation pressure and Poynting-Robertson drag).
 * @param rebx pointer to the rebx_extras instance
 * @param source_index Index in the particles array of the body that is the radiation source.
 * @param c Speed of light.
 */
struct rebx_params_radiation_forces* rebx_add_radiation_forces(struct rebx_extras* rebx, int source_index, double c);

/**
 * @brief Adds mass loss/growth to the simulation.
 * @param rebx pointer to the rebx_extras instance
 */
void rebx_add_modify_mass(struct rebx_extras* rebx);

/**
 * @brief Allows user to specify their own post timestep modifications. Behavior is identical to setting up post timestep modifications in REBOUND itself.
 * @param rebx pointer to the rebx_extras instance
 * @param custom_post_timestep_modifications Custom post-timestep modification function.
 * @param custom_params Custom parameters container.  Pass NULL if you don't want to use it.
 */
void rebx_add_custom_post_timestep_modification(struct rebx_extras* rebx, void (*custom_post_timestep_modification)(struct reb_simulation* const sim, struct rebx_effect* const custom_effect), void* custom_params);

/**
 * @brief Allows user to specify their own extra forces. Behaviour is identical to setting up extra forces  in REBOUND itself.
 * @param rebx pointer to the rebx_extras instance
 * @param custom_forces Custom forces function.
 * @param force_is_velocity_dependent Set to 1 if force is velocity dependent.
 * @param custom_params Custom parameters container.  Pass NULL if you don't want to use it.
 */
void rebx_add_custom_force(struct rebx_extras* rebx, void (*custom_force)(struct reb_simulation* const sim, struct rebx_effect* const custom_effect), int force_is_velocity_dependent, void* custom_params);
/** @} */
/** @} */

/********************************
 * Parameter getters and setters
 *******************************/

/**
 * \name Getters and Setters for particle parameters
 * @{
 */
/**
 * @defgroup GetterSetter
 * @brief Getters and setters for particle parameters (one for each variable type).
 * @{
 */

/**
 * @brief Sets a parameter of type double for a particular particle.
 * @param p Pointer to the particle in the simulation to which we want to add the parameter.
 * @param param_name Name of the parameter we want to set (see Effects page at http://reboundx.readthedocs.org)
 * @param value Value to which we want to set the parameter.
 */
void rebx_set_param_double(struct reb_particle* p, const char* param_name, double value);

/**
 * @brief Gets the parameter value of a particular particle.
 * @param p Pointer to the particle we want the parameter value for.
 * @param param_name Name of the parameter we want to get (see Effects page at http://reboundx.readthedocs.org)
 */
double rebx_get_param_double(struct reb_particle* p, const char* param_name);

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
 * @param rebx pointer to the rebx_extras instance.
 * @param params parameters structure returned when adding effect.
 * @param particle_radius radius of grain.
 * @param density density of particle.
 * @param Q_pr Radiation pressure coefficient (Burns et al. 1979).
 * @param L Luminosity of radiation source.
 */
double rebx_rad_calc_beta(struct rebx_extras* rebx, struct rebx_params_radiation_forces* params, double particle_radius, double density, double Q_pr, double L);
/**
 * @brief Calculates the particle radius from physical parameters and beta, the ratio of radiation to gravitational forces from the star.
 * @param rebx pointer to the rebx_extras instance.
 * @param params parameters structure returned when adding effect.
 * @param beta ratio of radiation force to gravitational force from the radiation source body.
 * @param density density of particle.
 * @param Q_pr Radiation pressure coefficient (Burns et al. 1979).
 * @param L Luminosity of radiation source.
 */
double rebx_rad_calc_particle_radius(struct rebx_extras* rebx, struct rebx_params_radiation_forces* params, double beta, double density, double Q_pr, double L);

/**
 * @brief Calculates the hamiltonian for gr_potential, including the classical Hamiltonian.
 * @param sim pointer to the REBOUND simulation
 * @param params parameters structure returned by add_gr_potential.
 */
double rebx_gr_potential_hamiltonian(const struct reb_simulation* const sim, const struct rebx_params_gr_potential* const params);

/**
 * @brief Calculates the hamiltonian for gr, including the classical Hamiltonian.
 * @param sim pointer to the REBOUND simulation
 * @param params parameters structure returned by add_gr.
 */
double rebx_gr_hamiltonian(const struct reb_simulation* const sim, const struct rebx_params_gr* const params);

/**
 * @brief Calculates the hamiltonian for gr_full, including the classical Hamiltonian.
 * @param sim pointer to the REBOUND simulation
 * @param params parameters structure returned by add_gr_full.
 */
double rebx_gr_full_hamiltonian(const struct reb_simulation* const sim, const struct rebx_params_gr_full* const params);

/** @} */
/** @} */

#endif
