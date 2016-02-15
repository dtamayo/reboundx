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

/****************************************
Internal functions
*****************************************/

/* Main routines called each timestep. */
void rebx_forces(struct reb_simulation* sim);               // Calls all the forces that have been added to the simulation.
void rebx_post_timestep_modifications(struct reb_simulation* sim);                  // Calls all the post-timestep modifications that have been added to the simulation.

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx); // Initializes all pointers and values.

/* Garbage collection routines. */
void rebx_free_params(struct rebx_extras* rebx);            // Steps through linked list to free all allocated parameters.
void rebx_free_effects(struct rebx_extras* rebx);        // Frees all effects in effects linked list 
void rebx_free_pointers(struct rebx_extras* rebx);          // Frees all the remaining pointers in rebx_extras.

/* Internal utility functions. */
void* rebx_get_param(const struct reb_particle* p, uint32_t param_type);  // returns rebx_param corresponding to the passed param in the passed particle.  If it doesn't exist, returns NULL.
    
void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param); // add a node for param in the rebx_params_to_be_freed linked list.

/* Internal parameter adders (need a different one for each type of particle parameter). */
void rebx_add_param_double(struct reb_particle* p, uint32_t param_type, double value);

/** @} */

/* Function for testing whether REBOUNDx can load librebound.so and call REBOUND functions. */
double install_test(void);

/****************************************
  New hash based functions
*****************************************/

uint32_t rebx_hash(const char* str);
void rebx_set_param_double_hash(struct reb_particle* p, uint32_t h, double value);
double rebx_get_param_double_hash(struct reb_particle* p, uint32_t h);

void rebx_add_force(struct rebx_extras* rebx, void* paramsPtr, const char* name, void (*functionPtr) (struct reb_simulation* sim, struct rebx_effect* effect));
void rebx_add_post_timestep_modification(struct rebx_extras* rebx, void* paramsPtr, const char* name, void (*functionPtr) (struct reb_simulation* sim, struct rebx_effect* effect));


#endif
