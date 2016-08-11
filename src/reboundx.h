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
#include "rebxtools.h"

extern const char* rebx_build_str;      ///< Date and time build string.
extern const char* rebx_version_str;    ///<Version string.

/******************************************
  Enums that might be shared across effects
*******************************************/

enum rebx_object_type{ 
    REBX_TYPE_EFFECT,
    REBX_TYPE_PARTICLE,
};

enum rebx_param_type{
    REBX_TYPE_DOUBLE,
    REBX_TYPE_INT,
};

/**
 * @brief Enumeration for different coordinate systems.
 */
enum REBX_COORDINATES{
    REBX_COORDINATES_JACOBI,                        ///< Jacobi coordinates.  Default for REBOUND/REBOUNDx.
    REBX_COORDINATES_BARYCENTRIC,                   ///< Coordinates referenced to the center of mass of the whole system.
    REBX_COORDINATES_PARTICLE,                      ///< Coordinates referenced to a particular particle.
};

/****************************************
Basic types in REBOUNDx
*****************************************/

/* 	Main structure used for all parameters added to particles.
 	These get added as nodes to a linked list for each particle, stored at particles[i].ap.*/
struct rebx_param{
    void* contents;                     // Pointer to the parameter (void* so it can point to different types).
    uint32_t hash;                      // Hash for the parameter name.
    enum rebx_param_type param_type;    // Enum for the parameter type.
    unsigned int length;                // 1 if paramPtr points to single value, more if an array.
    struct rebx_param* next;            // Pointer to the next parameter in the linked list.
};

/*  Structure for all REBOUNDx effects.
 *  These get added as nodes to the effects linked list in the rebx_extras structure.*/
struct rebx_effect{
    uint32_t hash;                      // hash corresponding to the effect's name.
    struct rebx_param* ap;              // Linked list of parameters for the effect.
    void (*force) (struct reb_simulation* sim, struct rebx_effect* effect); // Pointer to function to call during forces evaluation.
    void (*ptm) (struct reb_simulation* sim, struct rebx_effect* effect);   // Pointer to function to call after each timestep.
    struct rebx_extras* rebx;           // Pointer to the rebx_extras instance effect is in.
	struct rebx_effect* next;			// Pointer to the next effect in the linked list.
};

/*	Nodes for a linked list to all the parameters that have been allocated by REBOUNDx (so it can later free them).*/
struct rebx_param_to_be_freed{
    struct rebx_param* param;           // Pointer to a parameter node allocated by REBOUNDx.
    struct rebx_param_to_be_freed* next;// Pointer to the next node in the linked list rebx_extras.params_to_be_freed.
};

/*  Main REBOUNDx structure*/
struct rebx_extras {	
	struct reb_simulation* sim;								// Pointer to the simulation REBOUNDx is linked to.
	struct rebx_effect* effects;		                    // Linked list with pointers to all the effects added to the simulation.
	struct rebx_param_to_be_freed* params_to_be_freed; 		// Linked list with pointers to all parameters allocated by REBOUNDx (for later freeing).
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
 * @brief Removes a parameter from a particle.
 * @param effect Pointer to the particle we want to remove a parameter from.
 * @param param_name Name of the parameter we want to remove.
 * @return 1 if parameter found and successfully removed, 0 otherwise.
 */
int rebx_remove_param(const char* const param_name, const void* const object, enum rebx_object_type object_type);

/**
 * @brief Sets a parameter of type double for a particle.
 * @param object Pointer to the particle to which to add the parameter.
 * @param param_name Name of the parameter we want to set (see Effects page at http://reboundx.readthedocs.org for what parameters are needed for each effect)
 * @param value Value to which we want to set the parameter.
 */
void rebx_set_param(const char* const param_name, void* const value, enum rebx_param_type param_type, const unsigned int length, const void* const object, enum rebx_object_type object_type);

/**
 * @brief Gets a parameter value of type double from a REBOUNDx effect.
 * @param object Pointer to the effect that holds the parameter.
 * @param param_name Name of the parameter we want to get (see Effects page at http://reboundx.readthedocs.org)
 * @return Pointer to the parameter. NULL if parameter is not found in object (user must check for NULL to avoid segmentation fault).
 */
int rebx_get_param(const char* const param_name, void* value, enum rebx_param_type param_type, const unsigned int length, const void* const object, enum rebx_object_type object_type);

double rebx_get_doubleP(const char* const param_name, const void* const object);
double rebx_get_doubleE(const char* const param_name, const void* const object);
void rebx_set_doubleP(const char* const param_name, double value, const void* const object);
void rebx_set_doubleE(const char* const param_name, double value, const void* const object);
void rebx_set_doublesP(const char* const param_name, double* const value, const unsigned int length, const void* const object);
void rebx_get_doublesP(const char* const param_name, double* const value, const unsigned int length, const void* const object);
void rebx_set_doublesE(const char* const param_name, double* const value, const unsigned int length, const void* const object);
void rebx_get_doublesE(const char* const param_name, double* const value, const unsigned int length, const void* const object);
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
 * @brief Calculates the Aradial parameter for central_force effect required for a particle to have a particular pericenter precession rate.
 * @param p Particle whose pericenter precession rate we want to match.
 * @param primary Central particle for the central force (to which we add the Acentral and gammacentral parameters).
 * @param pomegadot Pericenter precession rate we want to obtain.
 * @param gamma Index of the central force law.
 * @return Acentral Normalization to add to the central particle.
 */
double rebx_central_force_Acentral(const struct reb_particle p, const struct reb_particle primary, const double pomegadot, const double gamma);

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

/**
 * @brief Calculates the hamiltonian for tides_precession effect.
 * @param sim pointer to the REBOUND simulation
 * @return Potential corresponding to tides_precession effect.
 */
double rebx_tides_precession_hamiltonian(struct reb_simulation* const sim);

/**
 * @brief Calculates the hamiltonian for central_force effect.
 * @param sim pointer to the REBOUND simulation.
 * @return Potential corresponding to central_force effect.
 */
double rebx_central_force_hamiltonian(struct reb_simulation* const sim);

/** @} */
/** @} */

#endif
