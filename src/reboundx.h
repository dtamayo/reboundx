/**
 * @file    reboundx.h
 * @brief   REBOUNDx API definition.
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>, Hanno Rein
 *
 * @section     LICENSE
 * Copyright (c) 2019 Dan Tamayo, Hanno Rein
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

#include <stdint.h>
#include <limits.h>
#include "rebound.h"
#include "rebxtools.h"
#ifndef REBXGITHASH
#define REBXGITHASH notavailable0000000000000000000000000001
#endif // REBXGITHASH

extern const char* rebx_build_str;      ///< Date and time build string.
extern const char* rebx_version_str;    ///< Version string.
extern const char* rebx_githash_str;    ///< Current git hash.

/******************************************
  REBOUNDx Enums
*******************************************/

/**
 * @brief Data types available in REBOUNDx for parameters.
 */
enum rebx_param_type{
    REBX_TYPE_NONE,
    REBX_TYPE_DOUBLE,
    REBX_TYPE_INT,
    REBX_TYPE_POINTER,
    REBX_TYPE_FORCE,
    REBX_TYPE_UINT32,
    REBX_TYPE_ORBIT,
    REBX_TYPE_ODE,
    REBX_TYPE_VEC3D
};

/**
 * @brief Different coordinate systems.
 */
enum REBX_COORDINATES{
    REBX_COORDINATES_JACOBI,                ///< Jacobi coordinates (default)
    REBX_COORDINATES_BARYCENTRIC,           ///< Coordinates referenced to pos/vel of system's center of mass.
    REBX_COORDINATES_PARTICLE,              ///< Coordinates relative to pos/vel of a particular particle.
};

/**
 * @brief Flag for whether steps should happen before or after the timestep
 */
enum rebx_timing {
    REBX_TIMING_PRE = -1,   ///< Pre timestep
    REBX_TIMING_POST = 1,   ///< Post timestep
};

/**
 * @brief Force types
 */
enum rebx_force_type{
    REBX_FORCE_NONE,        ///< Uninitialized default
    REBX_FORCE_POS,         ///< Force derivable from a position-dependent potential
    REBX_FORCE_VEL,         ///< velocity (or pos and vel) dependent force
};

/**
 * @brief Operator types
 */
enum rebx_operator_type{
    REBX_OPERATOR_NONE,     ///< Uninitialized default
    REBX_OPERATOR_UPDATER,  ///< operator that modifies x,v or m,
    REBX_OPERATOR_RECORDER, ///< operator that leaves state unchanged. Just records
};

/**
 * @brief Different fields for binary files
 */
enum rebx_binary_field_type{
    REBX_BINARY_FIELD_TYPE_NONE=0,
    REBX_BINARY_FIELD_TYPE_OPERATOR=1,
    REBX_BINARY_FIELD_TYPE_PARTICLE=2,
    REBX_BINARY_FIELD_TYPE_REBX_STRUCTURE=3,
    REBX_BINARY_FIELD_TYPE_PARAM=4,
    REBX_BINARY_FIELD_TYPE_NAME=5,
    REBX_BINARY_FIELD_TYPE_PARAM_TYPE=6,
    REBX_BINARY_FIELD_TYPE_PARAM_VALUE=7,
    REBX_BINARY_FIELD_TYPE_END=8,
    REBX_BINARY_FIELD_TYPE_PARTICLE_INDEX=9,
    REBX_BINARY_FIELD_TYPE_REBX_INTEGRATOR=10,
    REBX_BINARY_FIELD_TYPE_FORCE_TYPE=11,
    REBX_BINARY_FIELD_TYPE_OPERATOR_TYPE=12,
    REBX_BINARY_FIELD_TYPE_STEP=13,
    REBX_BINARY_FIELD_TYPE_STEP_DT_FRACTION=14,
    REBX_BINARY_FIELD_TYPE_REGISTERED_PARAM=15,
    REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCE=16,
    REBX_BINARY_FIELD_TYPE_PARAM_LIST=17,
    REBX_BINARY_FIELD_TYPE_REGISTERED_PARAMETERS=18,
    REBX_BINARY_FIELD_TYPE_ALLOCATED_FORCES=19,
    REBX_BINARY_FIELD_TYPE_ALLOCATED_OPERATORS=20,
    REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCES=21,
    REBX_BINARY_FIELD_TYPE_PRE_TIMESTEP_MODIFICATIONS=22,
    REBX_BINARY_FIELD_TYPE_POST_TIMESTEP_MODIFICATIONS=23,
    REBX_BINARY_FIELD_TYPE_PARTICLES=24,
    REBX_BINARY_FIELD_TYPE_FORCE=25,
    REBX_BINARY_FIELD_TYPE_SNAPSHOT=26,
};

/**
 * @brief Possible errors that might occur during binary file reading.
 */

// not loaded should be warning, since user might have saved with custom params or effects
enum rebx_input_binary_messages {
    REBX_INPUT_BINARY_WARNING_NONE = 0,
    REBX_INPUT_BINARY_ERROR_NOFILE = 1,
    REBX_INPUT_BINARY_ERROR_CORRUPT = 2,
    REBX_INPUT_BINARY_ERROR_NO_MEMORY = 4,
    REBX_INPUT_BINARY_ERROR_REBX_NOT_LOADED = 8,
    REBX_INPUT_BINARY_ERROR_REGISTERED_PARAM_NOT_LOADED = 16,
    REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED = 32,
    REBX_INPUT_BINARY_WARNING_PARTICLE_PARAMS_NOT_LOADED = 64,
    REBX_INPUT_BINARY_WARNING_FORCE_NOT_LOADED = 128,
    REBX_INPUT_BINARY_WARNING_OPERATOR_NOT_LOADED = 256,
    REBX_INPUT_BINARY_WARNING_STEP_NOT_LOADED = 512,
    REBX_INPUT_BINARY_WARNING_ADDITIONAL_FORCE_NOT_LOADED = 1024,
    REBX_INPUT_BINARY_WARNING_FIELD_UNKNOWN = 2048,
    REBX_INPUT_BINARY_WARNING_LIST_UNKNOWN = 4096,
    REBX_INPUT_BINARY_WARNING_PARAM_VALUE_NULL = 8192,
    REBX_INPUT_BINARY_WARNING_VERSION = 16384,
    REBX_INPUT_BINARY_WARNING_FORCE_PARAM_NOT_LOADED = 32768,
};

/**
 * @brief Different schemes for integrating across the interaction step
 */
enum rebx_integrator {
    REBX_INTEGRATOR_NONE = -1,
    REBX_INTEGRATOR_IMPLICIT_MIDPOINT = 0,
    REBX_INTEGRATOR_RK4 = 1,
    REBX_INTEGRATOR_EULER = 2,
    REBX_INTEGRATOR_RK2 = 3,
};

/**
 * @brief Different interpolation options
 */
enum rebx_interpolation_type {
    REBX_INTERPOLATION_NONE = 0,
    REBX_INTERPOLATION_SPLINE = 1,
};

/****************************************
Basic types in REBOUNDx
*****************************************/

/**
 * @brief Node structure for all REBOUNDx linked lists.
 */

struct rebx_node{
    void* object;             ///< Pointer to object (param, force, step, etc)
    struct rebx_node* next;   ///< Pointer to next node in list
};

/**
 * @brief Main structure used for all parameters added to objects.
 */

struct rebx_param{
    char* name;                 ///< For searching linked lists and informative errors
    enum rebx_param_type type;  ///< Needed to cast value
    void* value;                ///< Pointer to parameter value
};

/**
 * @brief Structure for REBOUNDx forces.
 */
struct rebx_force{
    char* name;                 ///< For searching linked lists and informative errors
    struct rebx_node* ap;       ///< Additional parameters linked list
    struct reb_simulation* sim; ///< Pointer to attached sim. Needed for error checks
    // See comments in params.py in __init__
    enum rebx_force_type force_type;    ///< Force type for internal logic
    void (*update_accelerations) (struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N); ///< Function pointer to add additional accelerations
};

/**
 * @brief Structure for REBOUNDx operators.
 */
struct rebx_operator{
    char* name;                 ///< For searching linked lists and informative errors
    struct rebx_node* ap;       ///< Additional parameters linked list
    struct reb_simulation* sim; ///< Pointer to attached sim. Needed for error checks
    // See comments in params.py in __init__
    enum rebx_operator_type operator_type;  ///< Operator type for internal logic
    void (*step_function) (struct reb_simulation* sim, struct rebx_operator* operator, const double dt);       ///< Function pointer to execute step
};

/**
 * @brief Structure for a REBOUNDx step.
 * @details A step is just a combination of an operator with a fraction of a timestep (see Sec. 6 of REBOUNdx paper). Can use same operator for different steps of different lengths to build higher order splitting schemes.
 */
struct rebx_step{
    struct rebx_operator* operator;     ///< Pointer to operator to use
    double dt_fraction;                 ///< Fraction of sim.dt to use each time it's called
};

/**
 * @brief Structure used as building block to save and load binary files.
 */
struct rebx_binary_field{
    enum rebx_binary_field_type type;   ///< Type of object
    long size;                          ///< Size in bytes of the object data (not including this structure). So you can skip ahead.
};

struct rebx_interpolator{
    enum rebx_interpolation_type interpolation;
    double* times;
    double* values;
    int Nvalues;
    double* y2;
    int klo;
};
/**
 * @brief Main REBOUNDx structure.
 * @details These fields are used internally by REBOUNDx and generally should not be changed manually by the user. Use the API instead.
 */
struct rebx_extras {
	struct reb_simulation* sim;					    ///< Pointer to the simulation REBOUNDx is linked to.

    struct rebx_node* additional_forces;            ///< Linked list of extra forces
    struct rebx_node* pre_timestep_modifications;   ///< Linked list of rebx_steps to apply before each timestep
	struct rebx_node* post_timestep_modifications;  ///< Linked list of rebx_steps to apply after each timestep

    struct rebx_node* registered_params;            ///< Linked list of rebx_params with all the parameter names registered with their type (for type safety)
    struct rebx_node* allocated_forces;             ///< For memory management
    struct rebx_node* allocated_operators;          ///< For memory management
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
 * @param sim Pointer to the reb_simulation on which to add REBOUNDx functionality.
 * @return Pointer to a rebx_extras structure.
 */
struct rebx_extras* rebx_attach(struct reb_simulation* sim);

/**
 * @brief Detaches REBOUNDx from simulation, resetting all the simulation's function pointers that REBOUNDx has set.
 * @details This does not free the memory allocated by REBOUNDx (call rebx_free).
 * @param sim Pointer to the simulation from which to remove REBOUNDx
 */
void rebx_detach(struct reb_simulation* sim, struct rebx_extras* rebx);
void rebx_extras_cleanup(struct reb_simulation* sim);
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
void rebx_init_extras_from_binary(struct rebx_extras* rebx, const char* const filename, enum rebx_input_binary_messages* warnings);
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
int rebx_add_operator_step(struct rebx_extras* rebx, struct rebx_operator* operator, const double dt_fraction, enum rebx_timing timing);
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
 * @brief Get a pointer to a force by name.
 * @param rebx Pointer to the rebx_extras instance
 * @param effect_name Name of the force (string)
 * @return Pointer to the corresponding rebx_force structure, or NULL if not found.
 */
struct rebx_force* rebx_get_force(struct rebx_extras* const rebx, const char* const name);

/**
 * @brief Get a pointer to an operator by name.
 * @param rebx Pointer to the rebx_extras instance
 * @param effect_name Name of the operator (string)
 * @return Pointer to the corresponding rebx_operator structure, or NULL if not found.
 */
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
 * @details These functions make up the main interface for users.  See below for more specialized parameter functions.
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
 * @brief Gets a parameter from a particle or effect.
 * @param ap Pointer from which to get the param
 * @param param_name Name of the parameter we want to get (see Effects page at http://reboundx.readthedocs.org)
 * @return A void pointer to the parameter. NULL if not found.
 */

void* rebx_get_param(struct rebx_extras* const rebx, struct rebx_node* ap, const char* const param_name);
struct rebx_param* rebx_get_param_struct(struct rebx_extras* const rebx, struct rebx_node* ap, const char* const param_name);
void rebx_set_param_pointer(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, void* val);
void rebx_set_param_double(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, double val);
void rebx_set_param_int(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, int val);
void rebx_set_param_uint32(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, uint32_t val);
void rebx_set_param_vec3d(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, struct reb_vec3d val);
void rebx_register_param(struct rebx_extras* const rebx, const char* name, enum rebx_param_type type);

/** @} */
/** @} */

/******************************************
  General utility functions
*******************************************/

/**
 * \name General utility functions
 * @{
 */
/**
 * @defgroup UtilFunc
 * @brief General utility functions
 * @{
 */

/**
 * @brief Calculate spin angular momentum in the simulation of any bodies with spin parameters set (moment of inertia I and angular rotation frequency vector Omega).
 *
 * @param rebx Pointer to the rebx_extras instance
 */
struct reb_vec3d rebx_tools_spin_angular_momentum(struct rebx_extras* const rebx);

void rebx_simulation_irotate(struct rebx_extras* const rebx, const struct reb_rotation q);



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
 * @brief Count how many particles have their moment of inertia and spin set, and initialize corresponding spin ODEs
 * @details Must be called after setting the moment of inertia and spin of all particles you want to evolve. Attaches spin_ode param to passed effect struct.
 * @param rebx Pointer to the rebx_extras instance.
 * @param effect (rebx_force) Force structure to which to attach spin_ode object.
 */
void rebx_spin_initialize_ode(struct rebx_extras* const rebx, struct rebx_force* const effect);

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
 * @brief Calculates the hamiltonian for gr, including the classical Hamiltonian.
 * @details Assumes there is only one source particle (with gr_source set to 1)
 * @param rebx pointer to the REBOUNDx extras instance.
 * @param gr Force structure returned by rebx_load_force
 * @return Total Hamiltonian, including classical Hamiltonian (double).
 */
double rebx_gr_hamiltonian(struct rebx_extras* const rebx, const struct rebx_force* const gr);

/**
 * @brief Calculates the hamiltonian for gr_full, including the classical Hamiltonian.
 * @param rebx pointer to the REBOUNDx extras instance.
 * @param gr_full Force structure returned by rebx_load_force
 * @return Total Hamiltonian, including classical Hamiltonian (double).
 */
double rebx_gr_full_hamiltonian(struct rebx_extras* const rebx, const struct rebx_force* const gr_full);

/**
 * @brief Calculates the potential for gr_potential, including the classical Hamiltonian.
 * @param rebx pointer to the REBOUNDx extras instance.
 * @param gr_potential Force structure returned by rebx_load_force
 * @return Total Hamiltonian, including classical Hamiltonian (double).
 */
double rebx_gr_potential_potential(struct rebx_extras* const rebx, const struct rebx_force* const gr_potential);

/**
 * @brief Calculates the potential for the conservative piece of the tides_constant_time_lag effect.
 * @details Will be conserved if tctl_tau = 0
 * @param rebx pointer to the REBOUNDx extras instance.
 * @return Potential corresponding to tides_constant_time_lag effect.
 */
double rebx_tides_constant_time_lag_potential(struct rebx_extras* const rebx);

/**
 * @brief Calculates the total energy associated with bodies' spin and their tidal interactions in tides_spin effect. This includes the spin kinetic energy plus gravitational potential between the tidally and rotationally induced quadrupoles and other bodies (treated as point masses). You should add sim.energy() to include bodies' overall kinetic energy and point-mass gravitational potential energy.
 * @param rebx pointer to the REBOUNDx extras instance.
 * @return Energy associated with tides_spin effect.
 */
double rebx_tides_spin_energy(struct rebx_extras* const rebx);

/**
 * @brief Calculates the potential for central_force effect.
 * @param rebx pointer to the REBOUNDx extras instance.
 * @return Potential corresponding to central_force effect.
 */
double rebx_central_force_potential(struct rebx_extras* const rebx);

/**
 * @brief Calculates the potential for all particles with additional gravity field harmonics beyond the monopole (i.e., J2, J4).
 * @param rebx pointer to the REBOUNDx extras instance.
 * @return Potential corresponding to the effect from all particles of their additional gravity field harmonics
 */
double rebx_gravitational_harmonics_potential(struct rebx_extras* const rebx);

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
 * @details The main parameter manipulation functions are listed above, but these functions can be useful when writing effects.
 * @{
 */

/**
 * @brief Gets the full rebx_param structure for a particular parameter, rather than just the pointer to the contents.
 * @details This can be useful to check properties of the parameter, like the param_type or shape.
 * @param object Pointer to the particle or effect that holds the parameter.
 * @param param_name Name of the parameter we want to get (see Effects page at http://reboundx.readthedocs.org)
 * @return Pointer to the rebx_param structure that holds the parameter. NULL if not found.
 */
/**
 * @brief Returns a void pointer to the parameter just like rebx_get_param, but additionally checks that the param_type matches what is expected.
 * @details Effects should use this function rather than rebx_get_param to ensure that the user appropriately set parameters if working from Python.
 * @param object Pointer to the particle or effect that holds the parameter.
 * @param param_name Name of the parameter we want to get (see Effects page at http://reboundx.readthedocs.org)
 * @return Void pointer to the parameter. NULL if not found or type does not match (will write error to stderr).
 */
void* rebx_get_param_check(struct reb_simulation* sim, struct rebx_node* ap, const char* const param_name, enum rebx_param_type param_type);
struct rebx_param* rebx_get_or_add_param(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name);


/****************************************
 Stepper Functions
 *****************************************/
/**
 * \name REBOUND Stepper Functions
 * @{
 */
/**
 * @defgroup StepperFunctions
 * @details Wrapper functions to execute a step using various REBOUND integrators. Can be used to make custom splitting schemes.
 * @{
 */

/**
 * @brief Executes a step for passed time dt using the IAS15 integrator in REBOUND.
 * @param sim Pointer to the simulation to step.
 * @param operator Unused pointer (kept for consistency with other operators). Can pass NULL.
 * @param dt timestep for which to step in simulation time units.
 */
void rebx_ias15_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
/**
 * @brief Executes a Kepler step for passed time dt using the WHFast integrator in REBOUND.
 * @details Will use the coordinates and other options set in sim.ri_whfast
 * @param sim Pointer to the simulation to step.
 * @param operator Unused pointer (kept for consistency with other operators). Can pass NULL.
 * @param dt timestep for which to step in simulation time units.
 */
void rebx_kepler_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
/**
 * @brief Executes a jump step for passed time dt using the WHFast integrator in REBOUND.
 * @details Will use the coordinates and other options set in sim.ri_whfast
 * @param sim Pointer to the simulation to step.
 * @param operator Unused pointer (kept for consistency with other operators). Can pass NULL.
 * @param dt timestep for which to step in simulation time units.
 */
void rebx_jump_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
/**
 * @brief Executes an interaction step for passed time dt using the WHFast integrator in REBOUND.
 * @details Will use the coordinates and other options set in sim.ri_whfast
 * @param sim Pointer to the simulation to step.
 * @param operator Unused pointer (kept for consistency with other operators). Can pass NULL.
 * @param dt timestep for which to step in simulation time units.
 */
void rebx_interaction_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
/**
 * @brief Executes a drift step for passed time dt using the leapfrog integrator in REBOUND.
 * @param sim Pointer to the simulation to step.
 * @param operator Unused pointer (kept for consistency with other operators). Can pass NULL.
 * @param dt timestep for which to step in simulation time units.
 */
void rebx_drift_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
/**
 * @brief Executes a kick step for passed time dt using the leapfrog integrator in REBOUND.
 * @param sim Pointer to the simulation to step.
 * @param operator Unused pointer (kept for consistency with other operators). Can pass NULL.
 * @param dt timestep for which to step in simulation time units.
 */
void rebx_kick_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
/** @} */
/** @} */

/****************************************
 Interpolation Routines
 *****************************************/
/**
 * \name Interpolation Routines
 * @{
 */
/**
 * @defgroup InterpolationFunctions
 * @details Functions for doing parameter interpolation.
 * @{
 */

/**
 * @brief Takes an array of times and corresponding array of values and returns a structure that allows interpolation of values at arbitrary times.
 * @details See the parameter interpolation examples in C and Python.
 * @param rebx pointer to the REBOUNDx extras instance.
 * @param Nvalues Length of times and values arrays (must be equal for both).
 * @param times Array of times at which the corresponding values are supplied.
 * @param values Array of values at each corresponding time.
 * @param interpolation Enum specifying the interpolation method. Defaults to spline.
 * @return Pointer to a rebx_interpolator structure. Call rebx_interpolate to get values.
 */
struct rebx_interpolator* rebx_create_interpolator(struct rebx_extras* const rebx, const int Nvalues, const double* times, const double* values, enum rebx_interpolation_type interpolation);
/**
 * @brief Frees the memory for a rebx_interpolator structure.
 */
void rebx_free_interpolator(struct rebx_interpolator* const interpolator);

/**
 * @brief Interpolate value at arbitrary times.
 * @details Need to first rebx_create_interpolator with an array of times and corresponding values to interpolate between. See parameter interpolation examples.
 * @param rebx Pointer to the REBOUNDx extras instance.
 * @param interpolator Pointer to the rebx_interpolator structure to interpolate from.
 * @param time Time at which to interpolate value.
 * @return Interpolated value at passed time.
 */
double rebx_interpolate(struct rebx_extras* const rebx, struct rebx_interpolator* const interpolator, const double time);
/** @} */
/** @} */

/****************************************
 Testing Functions
 *****************************************/
/**
 * \name Testing Functions
 * @{
 */
/**
 * @defgroup TestingFunctions
 * @details Functions for testing REBOUNDx
 * @{
 */

/**
 * @brief Loads a binary file, reads the header, and gives back the file pointer for manual reading
 * @param filename File to open
 * @param warnings Pointer to warnings enum to store warnings that come up
 * @return Returns a pointer to the binary file at the position following the header
 */
FILE* rebx_input_inspect_binary(const char* const filename, enum rebx_input_binary_messages* warnings);

/**
 * @brief Read the next field in a binary file
 * @param inf Pointer to the input file
 * @return Returns rebx_binary_field struct, initialized to 0 if read fails
 */
struct rebx_binary_field rebx_input_read_binary_field(FILE* inf);

/**
 * @brief Skip forward in binary file
 * @param inf Pointer to the input file
 * @param field_size Length by which to skip from current file position
 */
void rebx_input_skip_binary_field(FILE* inf, long field_size);

/** @} */
/** @} */

void rebx_error(struct rebx_extras* rebx, const char* const msg);
#endif
