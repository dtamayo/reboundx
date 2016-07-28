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
#include "spring.h"

const char* rebx_build_str = __DATE__ " " __TIME__; // Date and time build string. 
const char* rebx_version_str = "2.10.0";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.

/*****************************
 Initialization routines.
 ****************************/

struct rebx_extras* rebx_init(struct reb_simulation* sim){  // reboundx.h
    struct rebx_extras* rebx = malloc(sizeof(*rebx));
    rebx_initialize(sim, rebx);
    return rebx;
}

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx){
    sim->extras = rebx;
    rebx->sim = sim;

    if(sim->additional_forces || sim->post_timestep_modifications){
        reb_warning(sim, "sim->additional_forces or sim->post_timestep_modifications was already set.  If you want to use REBOUNDx, you should add custom effects through REBOUNDx also.  See http://reboundx.readthedocs.org/en/latest/c_examples.html#adding-custom-post-timestep-modifications-and-forces for a tutorial.");
    }

    sim->additional_forces = rebx_forces;
    sim->post_timestep_modifications = rebx_post_timestep_modifications;

	rebx->params_to_be_freed = NULL;
	rebx->effects = NULL;
}

/*****************************
 Garbage Collection Routines
 ****************************/

void rebx_remove_from_simulation(struct reb_simulation* sim){
    sim->additional_forces = NULL;
    sim->post_timestep_modifications = NULL;
}

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

static struct rebx_effect* rebx_add_effect(struct rebx_extras* rebx, const char* name){
    struct rebx_effect* effect = malloc(sizeof(*effect));
    
    effect->object_type = reb_hash("rebx_effect");
    effect->hash = reb_hash(name);
    effect->ap = NULL;
    effect->force = NULL;
    effect->ptm = NULL;
    effect->rebx = rebx;
    
    
    effect->next = rebx->effects;
    rebx->effects = effect;
    
    return effect;
}

struct rebx_effect* rebx_add(struct rebx_extras* rebx, const char* name){
    struct rebx_effect* effect = rebx_add_effect(rebx, name);
    struct reb_simulation* sim = rebx->sim;

    if(effect->hash == reb_hash("gr")){
        sim->force_is_velocity_dependent = 1;
        effect->force = rebx_gr;
    }
    else if (effect->hash == reb_hash("gr_full")){
        sim->force_is_velocity_dependent = 1;
        effect->force = rebx_gr_full;
    }
    else if (effect->hash == reb_hash("gr_potential")){
        effect->force = rebx_gr_potential;
    }
    else if (effect->hash == reb_hash("modify_mass")){
        effect->ptm = rebx_modify_mass;
    }
    else if (effect->hash == reb_hash("modify_orbits_direct")){
        effect->ptm = rebx_modify_orbits_direct;
    }
    else if (effect->hash == reb_hash("modify_orbits_forces")){
        sim->force_is_velocity_dependent = 1;
        effect->force = rebx_modify_orbits_forces;
    }
    else if (effect->hash == reb_hash("radiation_forces")){
        sim->force_is_velocity_dependent = 1;
        effect->force = rebx_radiation_forces;
    }
    else if (effect->hash == reb_hash("spring_forces")){
        effect->force = rebx_spring_forces;
    }
    else{
        char str[100]; 
        sprintf(str, "Effect '%s' passed to rebx_add not found.\n", name);
        reb_error(sim, str);
    }
    
    return effect;
}

struct rebx_effect* rebx_add_custom_force(struct rebx_extras* rebx, const char* name, void (*custom_force)(struct reb_simulation* const sim, struct rebx_effect* const effect), const int force_is_velocity_dependent){
    struct rebx_effect* effect = rebx_add_effect(rebx, name);
    effect->force = custom_force;
    struct reb_simulation* sim = rebx->sim;
    if(force_is_velocity_dependent){
        sim->force_is_velocity_dependent = 1;
    }
    return effect;
}

struct rebx_effect* rebx_add_custom_post_timestep_modification(struct rebx_extras* rebx, const char* name, void (*custom_ptm)(struct reb_simulation* const sim, struct rebx_effect* const effect)){
    struct rebx_effect* effect = rebx_add_effect(rebx, name);
    effect->ptm = custom_ptm;
    return effect;
}
    
void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param){
    struct rebx_param_to_be_freed* newparam = malloc(sizeof(*newparam));
    newparam->param = param;

    newparam->next = rebx->params_to_be_freed;
    rebx->params_to_be_freed = newparam;
}

/*********************************************************************************
 Internal getter/setter functions
 ********************************************************************************/

// Internal function that detects object type by casting void* to different pointer types.
static enum rebx_object_type rebx_get_object_type(const void* const object){
    struct rebx_effect* effect = (struct rebx_effect*)object;
    if (effect->object_type == reb_hash("rebx_effect")){
        return REBX_EFFECT;
    }
    else{
        return REBX_REB_PARTICLE;
    }
}

// Internal function that sets a new param in particle or effect linked list
static void* rebx_set_newparam(void* const object, struct rebx_param** const newparamptr){
    struct rebx_param* newparam = *newparamptr;
    enum rebx_object_type object_type = rebx_get_object_type(object);
    struct rebx_extras* rebx;
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
        default:
            fprintf(stderr, "Unsupported case in rebx_set_newparam.\n");
            exit(1);
    }
    
    rebx_add_param_to_be_freed(rebx, newparam);
    return newparam->paramPtr;
}

/*static struct reb_simulation* rebx_get_sim(const void* const object){
    enum rebx_object_type object_type = rebx_get_object_type(object);
    switch(object_type){
        case REBX_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            return effect->rebx->sim;
        }
        case REBX_REB_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            return p->sim;
        }
        default:
            fprintf(stderr, "Unsupported case in rebx_get_sim.\n");
            exit(1);
    }
} 

static void rebx_get_param_warning(const void* const object, const char* const param_name)
{
    struct reb_simulation* sim = rebx_get_sim(object);
    char str[200]; // TODO add check for whether using C
    sprintf(str, "Param %s was not found. Accessing pointer returned from rebx_get_param will cause segmentation fault (check for NULL).\n", param_name);
    reb_warning(sim,str); 
}
*/

/*********************************************************************************
 Type-independent functions for dealing with parameters
 ********************************************************************************/

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
        default:
            fprintf(stderr, "Unsupported case in rebx_get_param_hash.\n");
            exit(1);
    }    
    
    while(current != NULL){
        if(current->hash == hash){
            return current->paramPtr;
        }
        current = current->next;
    }
    return NULL;
}

struct rebx_param* rebx_get_paramptr_hash(const void* const object, uint32_t hash){
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
        default:
            fprintf(stderr, "Unsupported case in rebx_get_param_hash.\n");
            exit(1);
    }
    
    while(current != NULL){
        if(current->hash == hash){
            return current;
        }
        current = current->next;
    }
    return NULL;
}

int rebx_remove_param(const void* const object, const char* const param_name){
    uint32_t hash = reb_hash(param_name);
    struct rebx_param* current;
    enum rebx_object_type object_type = rebx_get_object_type(object);
    switch(object_type){
        case REBX_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            current = effect->ap;
            if(current->hash == hash){
                effect->ap = current->next;
                return 1;
            }
            break;
        }
        case REBX_REB_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            current = p->ap;
            if(current->hash == hash){
                p->ap = current->next;
                return 1;
            }
            break;
        }
        default:
            fprintf(stderr, "Unsupported case in rebx_remove_param.\n");
            exit(1);
    }    
  
    while(current->next != NULL){
        if(current->next->hash == hash){
            current->next = current->next->next;
            return 1;
        }
        current = current->next;
    }
    return 0;
}
            
/*********************************************************************************
 Getters and Setters for particle parameters (need new set for each variable type)
 ********************************************************************************/

double* rebx_get_param_double(const void* const object, const char* const param_name){
    uint32_t hash = reb_hash(param_name);
    void* voidptr = rebx_get_param_hash(object, hash);
    if (voidptr == NULL){
        return NULL;
    }
    else{
        return (double*)voidptr;
    }
}

void rebx_set_param_double(void* object, const char* const param_name, double value){
    uint32_t hash = reb_hash(param_name);
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
    newparam->type_hash = reb_hash("double");

    return (double*)rebx_set_newparam(object, &newparam);
}

int* rebx_get_param_int(const void* const object, const char* const param_name){
    uint32_t hash = reb_hash(param_name);
    void* voidptr = rebx_get_param_hash(object, hash);
    if (voidptr == NULL){
        return NULL;
    }
    else{
        return (int*)voidptr;
    }
}

void rebx_set_param_int(void* object, const char* const param_name, int value){
    uint32_t hash = reb_hash(param_name);
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
    newparam->type_hash = reb_hash("int");

    return (int*)rebx_set_newparam(object, &newparam);
}

struct rebx_spring* rebx_get_param_spring(const void* const object, const char* const param_name){
    uint32_t hash = reb_hash(param_name);
    void* voidptr = rebx_get_param_hash(object, hash);
    if (voidptr == NULL){
        return NULL;
    }
    else{
        return (struct rebx_spring*)voidptr;
    }
}

void rebx_set_param_spring(void* object, const char* const param_name, struct rebx_spring value){
    uint32_t hash = reb_hash(param_name);
    struct rebx_spring* ptr = rebx_get_param_hash(object, hash);
    if(ptr == NULL){
        ptr = rebx_add_param_spring(object, hash);  // add a new parameter if it doesn't exist
    }
    *ptr = value;                                   // update existing value
}

struct rebx_spring* rebx_add_param_spring(void* object, uint32_t hash){
    struct rebx_param* newparam = malloc(sizeof(*newparam));
    newparam->paramPtr = malloc(sizeof(struct rebx_spring));
    newparam->hash = hash;
    newparam->type_hash = reb_hash("spring");

    return (struct rebx_spring*)rebx_set_newparam(object, &newparam);
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
