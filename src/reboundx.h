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
#ifndef _REBX_REBOUNDX_H
#define _REBX_REBOUNDX_H

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
    REBX_TYPE_NONE,
    REBX_TYPE_DOUBLE,                               ///< C type double
    REBX_TYPE_INT,                                  ///< C type int
    REBX_TYPE_POINTER,
    REBX_TYPE_FORCE,
};
/*REBX_TYPE_UINT32,                               ///< C type uint32_t
    REBX_TYPE_ORBIT,                                ///< reb_orbit structure
    REBX_TYPE_LONGLONG,                             ///< To hold python 64 bit int
    REBX_TYPE_OPERATOR,
    REBX_TYPE_FORCE,
    
};*/

/**
 * @brief Enumeration for different coordinate systems.
 */
enum REBX_COORDINATES{
    REBX_COORDINATES_JACOBI,                        ///< Jacobi coordinates.  Default for REBOUND/REBOUNDx.
    REBX_COORDINATES_BARYCENTRIC,                   ///< Coordinates referenced to the center of mass of the whole system.
    REBX_COORDINATES_PARTICLE,                      ///< Coordinates referenced to a particular particle.
};

enum rebx_timing {
    REBX_TIMING_PRE = -1,
    REBX_TIMING_POST = 1,
};

enum rebx_force_type{
    REBX_FORCE_NONE,
    REBX_FORCE_POS,         // only position dependent force
    REBX_FORCE_VEL,         // velocity (or pos and vel) dependent force
};
enum rebx_operator_type{
    REBX_OPERATOR_NONE,
    REBX_OPERATOR_UPDATER,  // operator that modifies x,v or m,
    REBX_OPERATOR_RECORDER, // operator that leaves state unchanged. Just records
};


/****************************************
Basic types in REBOUNDx
*****************************************/

/**
 * @brief Node structure for all REBOUNDx linked lists.
 * @param name Name used for searching the linked list
 * @param object Data held by the node. Could be param, force, step etc
 * @param next Pointer to the next node in the linked list
 */

struct rebx_node{
    void* object;                   ///< Pointer to param
    struct rebx_node* next;   ///< Pointer to next node in list
};

/**
 * @brief Main structure used for all parameters added to particles.
 */

struct rebx_param{
    char* name;                     ///< Used to search linked lists and informative errors
    enum rebx_param_type type;      ///< Needed to cast value
    void* value;                    ///< Pointer to value
};

/**
 * @brief Structure for REBOUNDx operators.
 */
struct rebx_operator{
    char* name;
    struct rebx_node* ap;
    struct reb_simulation* sim;
    enum rebx_operator_type operator_type;
    void (*step) (struct reb_simulation* sim, struct rebx_operator* operator, const double dt);   ///< Pointer to function to call before and/or after each timestep.
};

/**
 * @brief Structure for REBOUNDx forces.
 */
struct rebx_force{
    char* name;
    struct rebx_node* ap;
    struct reb_simulation* sim; // need in python to match particles struct
    enum rebx_force_type force_type;
    void (*update_accelerations) (struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N); ///< Pointer to function to call during forces evaluation.
};

struct rebx_step{
    char* name;
    struct rebx_operator* operator;
    double dt_fraction;
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
    
    struct rebx_node* additional_forces;
    struct rebx_node* pre_timestep_modifications;		            ///< Linked list with pointers to all the effects added to the simulation.
	struct rebx_node* post_timestep_modifications;		            ///< Linked list with pointers to all the effects added to the simulation.
	
    struct rebx_node* registered_params;    ///< Linked list of all the names registered for parameters, along with their type
    struct rebx_node* allocated_forces;
    struct rebx_node* allocated_operators;
    
    enum {
        REBX_INTEGRATOR_IMPLICIT_MIDPOINT = 0,
        REBX_INTEGRATOR_RK4 = 1,
        REBX_INTEGRATOR_EULER = 2,
        REBX_INTEGRATOR_RK2 = 3,
        REBX_INTEGRATOR_NONE = -1,
    } integrator;
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
struct rebx_extras* rebx_attach(struct reb_simulation* sim);

/**
 * @brief Detaches REBOUNDx instance from simulation, resetting simulation's function pointers.
 * @param sim Pointer to the simulation from which to remove REBOUNDx
 */
void rebx_detach(struct reb_simulation* sim);

/**
 * @brief Frees all memory allocated by REBOUNDx instance.
 * @details Should be called after simulation is done if memory is a concern.
 * @param rebx The rebx_extras pointer returned from the initial call to rebx_attach.
 */
void rebx_free(struct rebx_extras* rebx);

int rebx_remove_force(struct rebx_extras* rebx, struct rebx_force* force);
int rebx_remove_operator(struct rebx_extras* rebx, struct rebx_operator* operator);

/**
 * @brief Save a binary file with all the effects in the simulation, as well as all particle and effect parameters.
 * @param rebx Pointer to the rebx_extras instance
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
 * @param rebx Pointer to the rebx_extras instance
 * @param name Name of the effect we want to add
 * @return Returns a pointer to a rebx_effect structure for the effect.
 */
//struct rebx_effect* rebx_add(struct rebx_extras* rebx, const char* name);
int rebx_add_operator(struct rebx_extras* rebx, struct rebx_operator* operator);
int rebx_add_operator_step(struct rebx_extras* rebx, struct rebx_operator* operator, const double dt_fraction, enum rebx_timing timing, char* name);
int rebx_add_force(struct rebx_extras* rebx, struct rebx_force* force);
struct rebx_operator* rebx_load_operator(struct rebx_extras* const rebx, const char* name);
struct rebx_force* rebx_load_force(struct rebx_extras* const rebx, const char* name);
struct rebx_operator* rebx_create_operator(struct rebx_extras* const rebx, const char* name);
struct rebx_force* rebx_create_force(struct rebx_extras* const rebx, const char* name);
/**
 * @brief Function for adding a custom force in REBOUNDx.
 * @param rebx Pointer to the rebx_extras instance
 * @param name String with the name of the custom effect
 * @param custom_force User-implemented function that updates the accelerations of particles.
 * @param force_is_velocity_dependent Should be set to 1 if the custom force uses particle velocities, 0 otherwise
 * @return Returns a pointer to a rebx_effect structure for the effect
 */
struct rebx_effect* rebx_add_custom_force(struct rebx_extras* rebx, const char* name, void (*custom_force)(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N), const int force_is_velocity_dependent);

/**
 * @brief Function for adding a custom post_timestep_modification in REBOUNDx.
 * @param rebx Pointer to the rebx_extras instance
 * @param name String with the name of the custom effect
 * @param custom_ptm User-implemented function that updates particles.
 * @return Returns a pointer to a rebx_effect structure for the effect.
 */
struct rebx_effect* rebx_add_custom_operator(struct rebx_extras* rebx, const char* name, void (*custom_operator)(struct reb_simulation* const sim, struct rebx_effect* const effect, const double dt, enum rebx_timing timing));
/**
 * @brief Get a pointer to an effect by name.
 * @param rebx Pointer to the rebx_extras instance
 * @param effect_name Name of the effect (string)
 * @return Pointer to the corresponding rebx_effect structure, or NULL if not found.
 */
struct rebx_force* rebx_get_force(struct rebx_extras* const rebx, const char* const name);
struct rebx_operator* rebx_get_operator(struct rebx_extras* const rebx, const char* const name);
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
int rebx_remove_param(struct rebx_node** apptr, const char* const param_name);

/**
 * @brief Adds a parameter to a particle or effect.
 * @param object Pointer to the particle or effect to which to add the parameter.
 * @param param_name Name of the parameter we want to set (see Effects page at http://reboundx.readthedocs.org for what parameters are needed for each effect)
 * @param param_type Variable type from rebx_param_type enumeration.
 * @return A void pointer to the parameter, i.e., the contents member of the new rebx_param structure. 
 */
//void* rebx_add_param(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, enum rebx_param_type param_type);
/**
 * @brief Adds an array parameter to a particle or effect.
 * @param object Pointer to the particle or effect to which to add the parameter.
 * @param param_name Name of the parameter we want to set (see Effects page at http://reboundx.readthedocs.org for what parameters are needed for each effect)
 * @param param_type Variable type from rebx_param_type enumeration.
 * @param ndim Number of dimensions in the array.
 * @param shape Pointer to an integer array specifying the length of the array in each dimension.
 * @return A void pointer to the parameter, i.e., the contents member of the new rebx_param structure.
 */
void* rebx_add_param_array(struct reb_simulation* const sim, struct rebx_node** apptr, const char* const param_name, enum rebx_param_type param_type, const int ndim, const int* const shape);


/**
 * @brief Gets a parameter from a particle or effect.
 * @param object Pointer to the particle or effect that holds the parameter.
 * @param param_name Name of the parameter we want to get (see Effects page at http://reboundx.readthedocs.org)
 * @return A void pointer to the parameter. NULL if not found.
 */

void* rebx_get_param(struct rebx_extras* rebx, struct rebx_node* ap, const char* const param_name);
struct rebx_param* rebx_get_param_struct(struct rebx_extras* rebx, struct rebx_node* ap, const char* const param_name);
int rebx_set_param_pointer(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, void* val);
int rebx_set_param_double(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, double val);
int rebx_set_param_int(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, int val);
void rebx_register_param(struct rebx_extras* const rebx, const char* name, enum rebx_param_type type);
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
double rebx_gr_hamiltonian(struct reb_simulation* const sim, const struct rebx_force* const gr);

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

/**
 * @brief Calculates the hamiltonian contribution for all particles with additional gravity field harmonics beyond the monopole (i.e., J2, J4).
 * @param sim pointer to the REBOUND simulation.
 * @return Potential corresponding to the effect from all particles of their additional gravity field harmonics
 */
double rebx_gravity_fields_hamiltonian(struct reb_simulation* const sim);

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
/**
 * @brief Returns a void pointer to the parameter just like rebx_get_param, but additionally checks that the param_type matches what is expected.
 * @detail Effects should use this function rather than rebx_get_param to ensure that the user appropriately set parameters if working from Python.
 * @param object Pointer to the particle or effect that holds the parameter.
 * @param param_name Name of the parameter we want to get (see Effects page at http://reboundx.readthedocs.org)
 * @return Void pointer to the parameter. NULL if not found or type does not match (will write error to stderr).
 */
void* rebx_get_param_check(struct reb_simulation* sim, struct rebx_node* ap, const char* const param_name, enum rebx_param_type param_type);

void rebx_gr_acc(struct rebx_extras* const rebx, double* acc, const double C2);
double rebx_calculate_energy(struct reb_simulation* const sim);
int rebx_len(struct rebx_node* head);
struct rebx_param_wrapper* rebx_get_param_wrapper(struct rebx_node* ap, const char* const param_name);

/*void rebx_ias15_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
void rebx_kepler_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
void rebx_jump_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
void rebx_interaction_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
*/
#endif
