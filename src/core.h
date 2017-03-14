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

// Nodes for a linked list to all the parameters that have been allocated by REBOUNDx (so it can later free them).
struct rebx_param_to_be_freed{
    struct rebx_param* param;           // Pointer to a parameter node allocated by REBOUNDx.
    struct rebx_param_to_be_freed* next;// Pointer to the next node in the linked list rebx_extras.params_to_be_freed.
};

/*****************************
 Internal initialization routine.
 ****************************/

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx); // Initializes all pointers and values.

/*****************************
 Garbage Collection Routines
 ****************************/

void rebx_free_params(struct rebx_extras* rebx);            // Steps through linked list to free all allocated particle parameters.
void rebx_free_effects(struct rebx_extras* rebx);           // Frees all effects in effects linked list 

/**********************************************
 Functions executing forces & ptm each timestep
 *********************************************/

void rebx_forces(struct reb_simulation* sim);                       // Calls all the forces that have been added to the simulation.
void rebx_pre_timestep_modifications(struct reb_simulation* sim);   // Calls all the pre-timestep modifications that have been added to the simulation.
void rebx_post_timestep_modifications(struct reb_simulation* sim);  // Calls all the post-timestep modifications that have been added to the simulation.

/**********************************************
 Adders 
 *********************************************/

// Add a parameter to the params_to_be_freed linked list for later freeing.
void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param); // add a node for param in the rebx_params_to_be_freed linked list.

/***********************************************************************************
 * Miscellaneous Functions
***********************************************************************************/
struct rebx_param* rebx_create_param();
struct rebx_param* rebx_attach_param_node(void* const object, struct rebx_param* param);
struct rebx_param* rebx_add_param_node(void* const object, const char* const param_name, enum rebx_param_type param_type, const int ndim, const int* const shape);
size_t rebx_sizeof(enum rebx_param_type param_type); // Returns size in bytes of the corresponding rebx_param_type type

double install_test(void);  // Function for testing whether REBOUNDx can load librebound.so and call REBOUND functions.

/****************************************
Force prototypes
*****************************************/
void rebx_gr(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N);
void rebx_gr_full(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N);
void rebx_gr_potential(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N);
void rebx_radiation_forces(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N);
void rebx_modify_orbits_forces(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N);
void rebx_tides_precession(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N);
void rebx_central_force(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N);
void rebx_tides_synchronous_ecc_damping(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N);

/****************************************
 Operator prototypes
 *****************************************/
void rebx_modify_mass(struct reb_simulation* const sim, struct rebx_effect* const effect, const double dt, enum rebx_timing timing);
void rebx_modify_orbits_direct(struct reb_simulation* const sim, struct rebx_effect* const effect, const double dt, enum rebx_timing timing);
void rebx_track_min_distance(struct reb_simulation* const sim, struct rebx_effect* const effect, const double dt, enum rebx_timing timing);

void rebx_reset_accelerations(struct reb_particle* const ps, const int N);
#endif
