/**
 * @file    core.h
 * @brief   Central internal functions for REBOUNDx (not called by user)
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

#ifndef _REBX_CORE_H
#define _REBX_CORE_H

#include <stdint.h>
#include "rebound.h"
#include "reboundx.h"

/**
 * @brief Enum for identifying different fields for binary files
 */
enum rebx_binary_field_type{
    REBX_BINARY_FIELD_TYPE_FORCE=0,
    REBX_BINARY_FIELD_TYPE_OPERATOR=1,
    REBX_BINARY_FIELD_TYPE_PARTICLE=2,
    REBX_BINARY_FIELD_TYPE_REBX_STRUCTURE=3,
    REBX_BINARY_FIELD_TYPE_PARAM=4,
    REBX_BINARY_FIELD_TYPE_NAME=5,
    REBX_BINARY_FIELD_TYPE_PARAM_TYPE=6,
    REBX_BINARY_FIELD_TYPE_VALUE=7,
    REBX_BINARY_FIELD_TYPE_END=8,
    REBX_BINARY_FIELD_TYPE_PARTICLE_INDEX=9,
    REBX_BINARY_FIELD_TYPE_REBX_INTEGRATOR=10,
    REBX_BINARY_FIELD_TYPE_FORCE_TYPE=11,
    REBX_BINARY_FIELD_TYPE_OPERATOR_TYPE=12,
    REBX_BINARY_FIELD_TYPE_STEP=13,
    REBX_BINARY_FIELD_TYPE_DT_FRACTION=14,
    REBX_BINARY_FIELD_TYPE_OPERATOR_NAME=15,
    REBX_BINARY_FIELD_TYPE_ADDITIONAL_FORCE=16,
};

/**
 * @brief Enum describing possible errors that might occur during binary file reading.
 */
enum rebx_input_binary_messages {
    REBX_INPUT_BINARY_WARNING_NONE = 0,
    REBX_INPUT_BINARY_ERROR_NOFILE = 1,
    REBX_INPUT_BINARY_ERROR_CORRUPT = 2,
    REBX_INPUT_BINARY_ERROR_NO_MEMORY = 4,
    REBX_INPUT_BINARY_WARNING_VERSION = 8,
    REBX_INPUT_BINARY_WARNING_PARAM_NOT_LOADED = 16,
    REBX_INPUT_BINARY_WARNING_PARTICLE_NOT_LOADED = 32,
    REBX_INPUT_BINARY_WARNING_FORCE_NOT_LOADED = 64,
    REBX_INPUT_BINARY_WARNING_OPERATOR_NOT_LOADED = 128,
    REBX_INPUT_BINARY_WARNING_STEP_NOT_LOADED = 256,
    REBX_INPUT_BINARY_WARNING_REG_PARAM_NOT_LOADED = 512,
    REBX_INPUT_BINARY_WARNING_FIELD_UNKOWN = 1024,
};

/*****************************
 Internal initialization routine.
 ****************************/

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx); // Initializes all pointers and values.

/**********************************************
 Functions executing forces & ptm each timestep
 *********************************************/

void rebx_additional_forces(struct reb_simulation* sim);                       // Calls all the forces that have been added to the simulation.
void rebx_pre_timestep_modifications(struct reb_simulation* sim);   // Calls all the pre-timestep modifications that have been added to the simulation.
void rebx_post_timestep_modifications(struct reb_simulation* sim);  // Calls all the post-timestep modifications that have been added to the simulation.

/**********************************************
 Adders 
 *********************************************/

// Add a parameter to the params_to_be_freed linked list for later freeing.
//void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param); // add a node for param in the rebx_params_to_be_freed linked list.

/***********************************************************************************
 * Miscellaneous Functions
***********************************************************************************/
//struct rebx_param* rebx_add_node(struct reb_simulation* const sim, struct rebx_param** head, const char* const param_name, enum rebx_param_type param_type, const int ndim, const int* const shape);
size_t rebx_sizeof(struct rebx_extras* rebx, enum rebx_param_type type); // Returns size in bytes of the corresponding rebx_param_type type
double install_test(void);  // Function for testing whether REBOUNDx can load librebound.so and call REBOUND functions.
void rebx_reset_accelerations(struct reb_particle* const ps, const int N);

/****************************************
Force prototypes
*****************************************/
void rebx_gr(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N);
void rebx_gr_full(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N);
void rebx_gr_potential(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N);
void rebx_radiation_forces(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N);
void rebx_modify_orbits_forces(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N);
void rebx_tides_precession(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N);
void rebx_central_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N);
void rebx_tides_synchronous_ecc_damping(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N);
void rebx_gravitational_harmonics(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N);

/****************************************
 Operator prototypes
 *****************************************/
void rebx_modify_mass(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
void rebx_modify_orbits_direct(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
void rebx_track_min_distance(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);

/****************************************
 Integrator prototypes
 *****************************************/

void rebx_integrator_euler_integrate(struct reb_simulation* const sim, const double dt, struct rebx_force* const force);
void rebx_integrator_rk2_integrate(struct reb_simulation* const sim, const double dt, struct rebx_effect* const effect);
void rebx_integrator_rk4_integrate(struct reb_simulation* const sim, const double dt, struct rebx_effect* const effect);
void rebx_integrator_implicit_midpoint_integrate(struct reb_simulation* const sim, const double dt, struct rebx_effect* const effect);

void* rebx_malloc(struct rebx_extras* const rebx, size_t memsize);
void rebx_free_ap(struct rebx_node** ap);
void rebx_free_particle_ap(struct reb_particle* p);
void rebx_free_force(struct rebx_force* force);
void rebx_free_operator(struct rebx_operator* operator);
void rebx_free_step(struct rebx_step* step);
void rebx_free_pointers(struct rebx_extras* rebx);
void rebx_free_param(struct rebx_param* param);

enum rebx_param_type rebx_get_type(struct rebx_extras* rebx, const char* name);

struct rebx_param* rebx_create_param(struct rebx_extras* rebx, const char* name, enum rebx_param_type type);
void rebx_add_param(struct rebx_extras* const rebx, struct rebx_node** apptr, struct rebx_param* param);
struct rebx_node* rebx_create_node(struct rebx_extras* rebx);

#endif
