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

#ifndef REBX_CORE_H
#define REBX_CORE_H

#include <stdint.h>
#include "rebound.h"

/****************************************
Basic types in REBOUNDx
*****************************************/

/* 	Main structure used for all parameters added to particles.
 	These get added as nodes to a linked list for each particle, stored at particles[i].ap.*/
struct rebx_param{
    void* paramPtr;                     // Pointer to the parameter (void* so it can point to different types of structs).
    uint32_t param_type;                // Identifier for the type of parameter.
    struct rebx_param* next;            // Pointer to the next parameter in the linked list.
};

/*  Structure for all REBOUNDx effects.
 *  These get added as nodes to the effects linked list in the rebx_extras structure.*/
struct rebx_effect{
	void* paramsPtr;					// Pointer to the effect params structure (void* so it can point to different effect structs).
	void (*functionPtr) (struct reb_simulation* const sim, struct rebx_effect* const effect);	// Pointer to the function to carry out the additional effect.
	uint32_t effect_type;       		// Identifier for the type of effect.
    int is_force;                       // Flag for whether effect is a force (1) or post_timestep_modification (0)
	struct rebx_effect* next;			// Pointer to the next effect in the linked list.
};

/*	Nodes for a linked list to all the parameters that have been allocated by REBOUNDx (so it can later free them).*/
struct rebx_param_to_be_freed{
    struct rebx_param* param;           // Pointer to a parameter node allocated by REBOUNDx.
    struct rebx_param_to_be_freed* next;// Pointer to the next node in the linked list rebx_extras.params_to_be_freed.
};

/****************************************
Main REBOUNDx structure
*****************************************/
struct rebx_extras {	
	struct reb_simulation* sim;								// Pointer to the simulation REBOUNDx is linked to.
	struct rebx_effect* effects;		                    // Linked list with pointers to all the effects added to the simulation.
	struct rebx_param_to_be_freed* params_to_be_freed; 		// Linked list with pointers to all parameters allocated by REBOUNDx (for later freeing).

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
void rebx_post_timestep_modifications(struct reb_simulation* sim);  // Calls all the post-timestep modifications that have been added to the simulation.

/**********************************************
 Adders for linked lists in extras structure
 *********************************************/

// Add a force to the effects linked list in the extras structure
void rebx_add_force(struct rebx_extras* rebx, void* paramsPtr, const char* name, void (*functionPtr) (struct reb_simulation* sim, struct rebx_effect* effect), int force_is_velocity_dependent);
// Add a post_timestep_modification to the effects linked list in the extras structure
void rebx_add_post_timestep_modification(struct rebx_extras* rebx, void* paramsPtr, const char* name, void (*functionPtr) (struct reb_simulation* sim, struct rebx_effect* effect));
// Add a parameter to the params_to_be_freed linked list for later freeing.
void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param); // add a node for param in the rebx_params_to_be_freed linked list.

/****************************************
 Hash functions
 *****************************************/

uint32_t rebx_murmur3_32(const char *key, uint32_t len, uint32_t seed); // hash function
uint32_t rebx_hash(const char* str);                                    // takes a string and returns a hash code

/*********************************************************************************
 General particle parameter getter
 ********************************************************************************/

void* rebx_get_param(const struct reb_particle* p, uint32_t param_type);// Generic parameter getter to be used by variable-type-specific functions below.  Returns rebx_param corresponding to the param_type hash code.  If it doesn't exist, returns NULL.

/*********************************************************************************
 Getters and Setters for particle parameters (need new set for each variable type)
 ********************************************************************************/

void rebx_set_param_double_hash(struct reb_particle* p, uint32_t h, double value);      // set parameter using hash code
void rebx_add_param_double(struct reb_particle* p, uint32_t param_type, double value);  // add a new parameter to linked list and set to value

double rebx_get_param_double_hash(struct reb_particle* p, uint32_t h);                  // get parameter using hash code

/****************************************
Custom Effect Adders
*****************************************/

// Add a custom post timestep modification.  See reboundx/examples/custom_ptm.
void rebx_add_custom_post_timestep_modification(struct rebx_extras* rebx, void (*custom_post_timestep_modification)(struct reb_simulation* const sim, struct rebx_effect* custom_effect), void* custom_params);

// Add a custom force.  See reboundx/examples/custom_ptm.
void rebx_add_custom_force(struct rebx_extras* rebx, void (*custom_force)(struct reb_simulation* const sim, struct rebx_effect* custom_effect), int force_is_velocity_dependent, void* custom_params);

/***********************************************************************************
 * Miscellaneous Functions
***********************************************************************************/

double install_test(void);  // Function for testing whether REBOUNDx can load librebound.so and call REBOUND functions. */

#endif
