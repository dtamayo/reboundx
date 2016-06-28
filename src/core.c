/**
 * @file    core.c
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

/* Main routines called each timestep. */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "core.h"
#include "rebound.h"
#include "gr.h"
#include "gr_potential.h"

const char* rebx_build_str = __DATE__ " " __TIME__; // Date and time build string. 
const char* rebx_version_str = "2.8.7";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.


/*****************************
 Initialization routines.
 ****************************/

struct rebx_extras* rebx_init(struct reb_simulation* sim){  // reboundx.h
    if(sim->additional_forces || sim->post_timestep_modifications){
        fprintf(stderr,"ERROR: sim->additional_forces or sim->post_timestep_modifications was already set.  If you want to use REBOUNDx, you need to add custom effects through REBOUNDx.  See http://reboundx.readthedocs.org/en/latest/c_examples.html#adding-custom-post-timestep-modifications-and-forces for a tutorial.");
        exit(1);
    }

    struct rebx_extras* rebx = malloc(sizeof(*rebx));
    rebx_initialize(sim, rebx);
    return rebx;
}

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx){
    sim->extras = rebx;
    rebx->sim = sim;

    sim->additional_forces = rebx_forces;
    sim->post_timestep_modifications = rebx_post_timestep_modifications;

	rebx->params_to_be_freed = NULL;
	rebx->effects = NULL;
}

/*****************************
 Garbage Collection Routines
 ****************************/

void rebx_free(struct rebx_extras* rebx){                   // reboundx.h
    rebx_free_params(rebx);
    rebx_free_effects(rebx);
    free(rebx);
}

void rebx_free_params(struct rebx_extras* rebx){
    struct rebx_param_to_be_freed* current = rebx->params_to_be_freed;
    struct rebx_param_to_be_freed* temp_next;
    while(current != NULL){
        temp_next = current->next;
        free(current->param->paramPtr);
        free(current->param);
        free(current);
        current = temp_next;
    }
}

void rebx_free_effects(struct rebx_extras* rebx){
    struct rebx_effect* current = rebx->effects;
    struct rebx_effect* temp_next;

    while(current != NULL){
        temp_next = current->next;
        free(current);
        current = temp_next;
    }
}

/**********************************************
 Functions executing forces & ptm each timestep
 *********************************************/

void rebx_forces(struct reb_simulation* sim){
	struct rebx_extras* rebx = sim->extras;
	struct rebx_effect* current = rebx->effects;
	while(current != NULL){
        if(current->force != NULL){
		    current->force(sim, current);
        }
		current = current->next;
	}
}

void rebx_post_timestep_modifications(struct reb_simulation* sim){
	struct rebx_extras* rebx = sim->extras;
	struct rebx_effect* current = rebx->effects;
	while(current != NULL){
        if(current->ptm != NULL){
		    current->ptm(sim, current);
        }
		current = current->next;
	}
}

/**********************************************
 Adders for linked lists in extras structure
 *********************************************/

struct rebx_effect* rebx_add_effect(struct rebx_extras* rebx, const char* name){
    struct rebx_effect* effect = malloc(sizeof(*effect));
    
    effect->object_type = reb_tools_hash("rebx_effect");
    effect->ap = NULL;
    effect->force = NULL;
    effect->ptm = NULL;
    effect->rebx = rebx;
    
    struct reb_simulation* sim = rebx->sim;
    uint32_t hash = reb_tools_hash(name);

    if(hash == reb_tools_hash("gr")){
        sim->force_is_velocity_dependent = 1;
        effect->force = rebx_gr;
    }
    /*else if (hash == reb_tools_hash("gr_potential")){
        effect->force = rebx_gr_potential;
    }*/
    else{
        fprintf(stderr, "Effect not found.\n");
        exit(1);
    }
    
    effect->next = rebx->effects;
    rebx->effects = effect;
    
    return effect;
}
    
void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param){
    struct rebx_param_to_be_freed* newparam = malloc(sizeof(*newparam));
    newparam->param = param;

    newparam->next = rebx->params_to_be_freed;
    rebx->params_to_be_freed = newparam;
}

/*********************************************************************************
 General particle parameter getter
 ********************************************************************************/
// Internal function that detects object type by casting void* to different pointer types.
static enum rebx_object_type rebx_get_object_type(const void* const object){
    struct rebx_effect* effect = (struct rebx_effect*)object;
    if (effect->object_type == reb_tools_hash("rebx_effect")){
        return REBX_EFFECT;
    }
    else{
        return REBX_REB_PARTICLE;
    }
}

void* rebx_get_param_hash(const void* const object, uint32_t hash){
    struct rebx_param* current;
    enum rebx_object_type object_type = rebx_get_object_type(object);
    switch(object_type){
        case REBX_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            current = effect->ap;
            break;
        }
        case REBX_REB_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            current = p->ap;
            break;
        }
    }    
    
    while(current != NULL){
        if(current->hash == hash){
            return current->paramPtr;
        }
        current = current->next;
    }
    return NULL;
}

int rebx_remove_param(const void* const object, const char* const param_name){
    uint32_t hash = reb_tools_hash(param_name);
    struct rebx_param* current;
    struct rebx_param** prev;
    enum rebx_object_type object_type = rebx_get_object_type(object);
    switch(object_type){
        case REBX_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            current = effect->ap;
            prev = &effect->ap;
            break;
        }
        case REBX_REB_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            current = p->ap;
            prev = &p->ap;
            break;
        }
    }    
   
    while(current != NULL){
        if(current->hash == hash){
            *prev = current->next;
            return 1;
        }
        prev = &current->next;
        current = current->next;
    }
    return 0;
}
            
/*********************************************************************************
 Getters and Setters for particle parameters (need new set for each variable type)
 ********************************************************************************/

int rebx_get_param_double(const void* const object, const char* const param_name, double* ptr){
    uint32_t hash = reb_tools_hash(param_name);
    void* voidptr = rebx_get_param_hash(object, hash);
    if (voidptr == NULL){
        return 0;
    }
    else{
        *ptr = *(double *)voidptr;
        return 1;
    }
}

void rebx_set_param_double(void* object, const char* const param_name, double value){
    uint32_t hash = reb_tools_hash(param_name);
    double* ptr = rebx_get_param_hash(object, hash);
    if(ptr == NULL){
        ptr = rebx_add_param_double(object, hash);  // add a new parameter if it doesn't exist
    }
    *ptr = value;                                   // update existing value
}

double* rebx_add_param_double(void* object, uint32_t hash){
    struct rebx_param* newparam = malloc(sizeof(*newparam));
    newparam->paramPtr = malloc(sizeof(double));
    newparam->hash = hash;
    newparam->type_hash = reb_tools_hash("double");

    struct rebx_extras* rebx;
    enum rebx_object_type object_type = rebx_get_object_type(object);
    switch(object_type){
        case REBX_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            newparam->next = effect->ap;
            effect->ap = newparam;
            rebx = effect->rebx;
            break;
        }
        case REBX_REB_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            newparam->next = p->ap;
            p->ap = newparam;
            rebx = p->sim->extras;
            break;
        }
    }
    
    rebx_add_param_to_be_freed(rebx, newparam);
    return newparam->paramPtr;
}

int rebx_get_param_int(const void* const object, const char* const param_name, int* ptr){
    uint32_t hash = reb_tools_hash(param_name);
    void* voidptr = rebx_get_param_hash(object, hash);
    if (voidptr == NULL){
        return 0;
    }
    else{
        *ptr = *(int *)voidptr;
        return 1;
    }
}

void rebx_set_param_int(void* object, const char* const param_name, int value){
    uint32_t hash = reb_tools_hash(param_name);
    int* ptr = rebx_get_param_hash(object, hash);
    if(ptr == NULL){
        ptr = rebx_add_param_int(object, hash);  // add a new parameter if it doesn't exist
    }
    *ptr = value;                                   // update existing value
}

int* rebx_add_param_int(void* object, uint32_t hash){
    struct rebx_param* newparam = malloc(sizeof(*newparam));
    newparam->paramPtr = malloc(sizeof(int));
    newparam->hash = hash;
    newparam->type_hash = reb_tools_hash("int");

    struct rebx_extras* rebx;
    enum rebx_object_type object_type = rebx_get_object_type(object);
    switch(object_type){
        case REBX_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            newparam->next = effect->ap;
            effect->ap = newparam;
            rebx = effect->rebx;
            break;
        }
        case REBX_REB_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            newparam->next = p->ap;
            p->ap = newparam;
            rebx = p->sim->extras;
            break;
        }
    }
    
    rebx_add_param_to_be_freed(rebx, newparam);
    return newparam->paramPtr;
}

/***********************************************************************************
 * Miscellaneous Functions
***********************************************************************************/

double install_test(void){
    struct reb_simulation* sim = reb_create_simulation();
    struct reb_particle p = {0};
    p.m = 1.; 
    reb_add(sim, p); 
    struct reb_particle p1 = reb_tools_orbit2d_to_particle(sim->G, p, 0., 1., 0.2, 0., 0.);
    reb_add(sim, p1);
    reb_integrate(sim, 1.);
    return sim->particles[1].x;
}
