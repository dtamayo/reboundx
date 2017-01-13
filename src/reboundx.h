/**
 * @file    reboundx.h
 * @brief   REBOUNDx API definition.
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>, Hanno Rein
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
#include <limits.h>
#include "rebound.h"
#include "rebxtools.h"
#ifndef REBXGITHASH
#define REBXGITHASH notavailable0000000000000000000000000001 
#endif // REBXGITHASH

extern const char* rebx_build_str;      ///< Date and time build string.
extern const char* rebx_version_str;    ///<Version string.
extern const char* rebx_githash_str;    ///< Current git hash.

/******************************************
  REBOUNDx Enums 
*******************************************/

/**
 * @brief Enumeration for the different data types available in REBOUNDx for parameters.
 */
enum rebx_param_type{
    REBX_TYPE_DOUBLE,                               ///< C type double
    REBX_TYPE_INT,                                  ///< C type int
    REBX_TYPE_UINT32,                               ///< C type uint32_t
    REBX_TYPE_ORBIT,                                ///< reb_orbit structure
    REBX_TYPE_LONGLONG,                             ///< To hold python 64 bit int
};

/**
 * @brief Enumeration for different coordinate systems.
 */
enum REBX_COORDINATES{
    REBX_COORDINATES_JACOBI,                        ///< Jacobi coordinates.  Default for REBOUND/REBOUNDx.
    REBX_COORDINATES_BARYCENTRIC,                   ///< Coordinates referenced to the center of mass of the whole system.
    REBX_COORDINATES_PARTICLE,                      ///< Coordinates referenced to a particular particle.
};

/**
 * @brief Internally used enum for identifying structs that can take parameters.
 */
enum rebx_object_type{
    REBX_OBJECT_TYPE_EFFECT=INT_MAX-2,
    REBX_OBJECT_TYPE_PARTICLE=INT_MAX-1,
};

/**
 * @brief Enum for identifying different fields for binary files
 */
enum rebx_binary_field_type{
    REBX_BINARY_FIELD_TYPE_EFFECT=0,
    REBX_BINARY_FIELD_TYPE_PARTICLE=1,
    REBX_BINARY_FIELD_TYPE_PARAM=2,
    REBX_BINARY_FIELD_TYPE_NAMELENGTH=3,
    REBX_BINARY_FIELD_TYPE_NAME=4,
    REBX_BINARY_FIELD_TYPE_PARAM_TYPE=5,
    REBX_BINARY_FIELD_TYPE_NDIM=6,
    REBX_BINARY_FIELD_TYPE_SHAPE=7,
    REBX_BINARY_FIELD_TYPE_CONTENTS=8,
    REBX_BINARY_FIELD_TYPE_END=9,
    REBX_BINARY_FIELD_TYPE_PARTICLE_INDEX=10,
    REBX_BINARY_FIELD_TYPE_PYTHON_TYPE=11,
};

/**
 * @brief Enum describing possible errors that might occur during binary file reading.
 */
enum rebx_input_binary_messages {
    REBX_INPUT_BINARY_WARNING_NONE = 0,
    REBX_INPUT_BINARY_ERROR_NOFILE = 1,
    REBX_INPUT_BINARY_ERROR_CORRUPT = 2,
    REBX_INPUT_BINARY_WARNING_VERSION = 4,
    REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED = 8,
    REBX_INPUT_BINARY_WARNING_PARTICLE_NOT_LOADED = 16,
    REBX_INPUT_BINARY_WARNING_EFFECT_NOT_LOADED = 32,
    REBX_INPUT_BINARY_WARNING_FIELD_UNKOWN = 64,
};


/****************************************
Basic types in REBOUNDx
*****************************************/

/**
 * @brief Main structure used for all parameters added to particles.
 */
struct rebx_param{
    void* contents;                                 ///< Pointer to the parameter (void* so it can point to different types).
    uint32_t hash;                                  ///< Hash for the parameter name.
    char* name;                                     ///< String for the parameter's name.
    enum rebx_param_type param_type;                ///< Enum for the parameter type.
    int python_type;                                ///< Used by python side to store python type
	int ndim;                                       ///< Number of dimensions (to support array parameters)
	int* shape;                                     ///< Array of length ndim for the array shape (NULL for scalars).
    int* strides;                                   ///< Strides along the different dimensions for array indexing (NULL for scalars).
	int size;                                       ///< Total number of values in array (1 for scalars).
    struct rebx_param* next;                        ///< Pointer to the next parameter in the linked list.
};

/**
 * @brief Structure for all REBOUNDx effects.
 * @detail These get added as nodes to the effects linked list in the rebx_extras structure.
 */
struct rebx_effect{
    uint32_t hash;                                  ///< Hash for the effect's name.
    char* name;                                     ///< String for the effect's name.
    struct rebx_param* ap;                          ///< Linked list of parameters for the effect.
    void (*force) (struct reb_simulation* sim, struct rebx_effect* effect); ///< Pointer to function to call during forces evaluation.
    void (*ptm) (struct reb_simulation* sim, struct rebx_effect* effect);   ///< Pointer to function to call after each timestep.
    struct rebx_extras* rebx;                       ///< Pointer to the rebx_extras instance effect is in.
	struct rebx_effect* next;			            ///< Pointer to the next effect in the linked list.
    char pad[100];                                  ///< Pad to be able to cast to reb_particle for get_object_type function.
};

/**
 * @brief This structure is used to save and load binary files.
 */
struct rebx_binary_field{
    enum rebx_binary_field_type type;               ///< Type of object
    long size;                                      ///< Size in bytes of the object data (not including this structure). So you can skip ahead.
};

/**
 * @brief Main REBOUNDx structure.
 */
struct rebx_extras {	
	struct reb_simulation* sim;					    ///< Pointer to the simulation REBOUNDx is linked to.
	struct rebx_effect* effects;		            ///< Linked list with pointers to all the effects added to the simulation.
	struct rebx_param_to_be_freed* params_to_be_freed; 	///< Linked list with pointers to all parameters allocated by REBOUNDx (for later freeing).
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
 * @param rebx The rebx_extras pointer returned from the initial call to rebx_init.
 */
void rebx_free(struct rebx_extras* rebx);

/**
 * @brief Save a binary file with all the effects in the simulation, as well as all particle and effect parameters.
 * @param rebx The rebx_extras pointer returned from the initial call to rebx_init.
 * @param filename Filename to which to save the binary file.
 */
void rebx_output_binary(struct rebx_extras* rebx, char* filename);

/**
 * @brief Reads a REBOUNDx binary file, loads all effects and parameters.
 * @param sim Pointer to the simulation to which the effects and parameters should be added.
 * @param filename Filename of the saved binary file.
 */
struct rebx_extras* rebx_create_extras_from_binary(struct reb_simulation* sim, const char* const filename);

/**
 * @brief Similar to rebx_create_extras_from_binary(), but takes an extras instance (must be attached to a simulation) and allows for manual message handling.
 * @param rebx Pointer to a rebx_extras instance to be updated.
 * @param filename Filename of the saved binary file.
 * @param warnings Pointer to an array of warnings to be populated during loading. 
 */
void rebx_create_extras_from_binary_with_messages(struct rebx_extras* rebx, const char* const filename, enum rebx_input_binary_messages* warnings);

/** @} */
/** @} */

/****************************************
  Functions for manipulating effects
*****************************************/
/**
 * \name Functions for manipulating effects in REBOUNDx
 * @{
 */
/**
 * @defgroup EffectManipulators
 * @details These are the functions for manipulating effects in REBOUNDx.
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

/**
 * @brief Get a pointer to an effect by name.
 * @param rebx Pointer to the rebx_extras instance returned by rebx_init.
 * @param effect_name Name of the effect (string).
 * @return Pointer to the corresponding rebx_effect structure, or NULL if not found.
 */
struct rebx_effect* rebx_get_effect(struct rebx_extras* const rebx, const char* const effect_name);

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
 * @detail These functions make up the main interface for users.  See below for more specialized parameter functions.
 * @{
 */

/**
 * @brief Removes a parameter from a particle or effect.
 * @param object Pointer to the particle or effect we want to remove a parameter from.
 * @param param_name Name of the parameter we want to remove.
 * @return 1 if parameter found and successfully removed, 0 otherwise.
 */
int rebx_remove_param(const void* const object, const char* const param_name);

/**
 * @brief Adds a parameter to a particle or effect.
 * @param object Pointer to the particle or effect to which to add the parameter.
 * @param param_name Name of the parameter we want to set (see Effects page at http://reboundx.readthedocs.org for what parameters are needed for each effect)
 * @param param_type Variable type from rebx_param_type enumeration.
 * @return A void pointer to the parameter, i.e., the contents member of the new rebx_param structure. 
 */
void* rebx_add_param(void* const object, const char* const param_name, enum rebx_param_type param_type);

/**
 * @brief Adds an array parameter to a particle or effect.
 * @param object Pointer to the particle or effect to which to add the parameter.
 * @param param_name Name of the parameter we want to set (see Effects page at http://reboundx.readthedocs.org for what parameters are needed for each effect)
 * @param param_type Variable type from rebx_param_type enumeration.
 * @param ndim Number of dimensions in the array.
 * @param shape Pointer to an integer array specifying the length of the array in each dimension.
 * @return A void pointer to the parameter, i.e., the contents member of the new rebx_param structure.
 */void* rebx_add_param_array(void* const object, const char* const param_name, enum rebx_param_type param_type, const int ndim, const int* const shape);


/**
 * @brief Gets a parameter from a particle or effect.
 * @param object Pointer to the particle or effect that holds the parameter.
 * @param param_name Name of the parameter we want to get (see Effects page at http://reboundx.readthedocs.org)
 * @return A void pointer to the parameter. NULL if not found.
 */
void* rebx_get_param(const void* const object, const char* const param_name);

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

double rebx_moon_quadrupole_laskar_hamiltonian(struct reb_simulation* const sim);

/** @} */
/** @} */

/********************************
 * Specialized Parameter manipulation functions
 *******************************/

/**
 * \name Specialized functions for accessing and modifying particle and effect parameters.
 * @{
 */
/**
 * @defgroup SpecializedParameterManipulators
 * @brief Specialized functions for accessing and modifying particle and effect parameters.
 * @detail The main parameter manipulation functions are listed above, but these functions can be useful when writing effects.
 * @{
 */

/**
 * @brief Gets the full rebx_param structure for a particular parameter, rather than just the pointer to the contents.
 * @detail This can be useful to check properties of the parameter, like the param_type or shape.
 * @param object Pointer to the particle or effect that holds the parameter.
 * @param param_name Name of the parameter we want to get (see Effects page at http://reboundx.readthedocs.org)
 * @return Pointer to the rebx_param structure that holds the parameter. NULL if not found.
 */
struct rebx_param* rebx_get_param_node(const void* const object, const char* const param_name);

/**
 * @brief Returns a void pointer to the parameter just like rebx_get_param, but additionally checks that the param_type matches what is expected.
 * @detail Effects should use this function rather than rebx_get_param to ensure that the user appropriately set parameters if working from Python.
 * @param object Pointer to the particle or effect that holds the parameter.
 * @param param_name Name of the parameter we want to get (see Effects page at http://reboundx.readthedocs.org)
 * @return Void pointer to the parameter. NULL if not found or type does not match (will write error to stderr).
 */
void* rebx_get_param_check(const void* const object, const char* const param_name, enum rebx_param_type param_type);

#endif
