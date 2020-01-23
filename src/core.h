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

struct rebx_extras;
struct rebx_param;
enum rebx_param_type;
struct rebx_step;
struct rebx_node;

#include <stdint.h>
#include "rebound.h"
#include "reboundx.h"



/*****************************
 Internal initialization routine.
 ****************************/

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx); // Initializes all pointers and values.
void rebx_register_default_params(struct rebx_extras* rebx); // Registers default params

/**********************************************
 Functions executing forces & ptm each timestep
 *********************************************/

void rebx_additional_forces(struct reb_simulation* sim);                       // Calls all the forces that have been added to the simulation.
void rebx_pre_timestep_modifications(struct reb_simulation* sim);   // Calls all the pre-timestep modifications that have been added to the simulation.
void rebx_post_timestep_modifications(struct reb_simulation* sim);  // Calls all the post-timestep modifications that have been added to the simulation.

/***********************************************************************************
 * Miscellaneous Functions
***********************************************************************************/
//struct rebx_param* rebx_add_node(struct reb_simulation* const sim, struct rebx_param** head, const char* const param_name, enum rebx_param_type param_type, const int ndim, const int* const shape);
size_t rebx_sizeof(struct rebx_extras* rebx, enum rebx_param_type type); // Returns size in bytes of the corresponding rebx_param_type type
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
void rebx_gravitational_harmonics(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N);
void rebx_ephemeris_forces(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N);

/****************************************
 Operator prototypes
 *****************************************/
void rebx_modify_mass(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
void rebx_integrate_force(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
void rebx_modify_orbits_direct(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);
void rebx_track_min_distance(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt);

/****************************************
 Integrator prototypes
 *****************************************/

void rebx_integrator_euler_integrate(struct reb_simulation* const sim, const double dt, struct rebx_force* const force);
void rebx_integrator_rk2_integrate(struct reb_simulation* const sim, const double dt, struct rebx_force* const force);
void rebx_integrator_rk4_integrate(struct reb_simulation* const sim, const double dt, struct rebx_force* const force);
void rebx_integrator_implicit_midpoint_integrate(struct reb_simulation* const sim, const double dt, struct rebx_force* const force);

void* rebx_malloc(struct rebx_extras* const rebx, size_t memsize);
void rebx_free_ap(struct rebx_node** ap);
void rebx_free_particle_ap(struct reb_particle* p);
void rebx_free_force(struct rebx_extras* rebx, struct rebx_force* force);
void rebx_free_operator(struct rebx_operator* operator);
void rebx_free_step(struct rebx_step* step);
void rebx_free_pointers(struct rebx_extras* rebx);
void rebx_free_param(struct rebx_param* param);

enum rebx_param_type rebx_get_type(struct rebx_extras* rebx, const char* name);

struct rebx_param* rebx_create_param(struct rebx_extras* rebx, const char* name, enum rebx_param_type type);
int rebx_add_param(struct rebx_extras* const rebx, struct rebx_node** apptr, struct rebx_param* param);
struct rebx_node* rebx_create_node(struct rebx_extras* rebx);

#endif
