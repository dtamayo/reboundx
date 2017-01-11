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
#define STRINGIFY(s) str(s)
#define str(s) #s

const char* rebx_build_str = __DATE__ " " __TIME__; // Date and time build string. 
const char* rebx_version_str = "2.16.0";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* rebx_githash_str = STRINGIFY(GITHASH);             // This line gets updated automatically. Do not edit manually.

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
	rebx->params_to_be_freed = NULL;
	rebx->effects = NULL;
    
    if(sim->additional_forces || sim->post_timestep_modifications){
        reb_warning(sim, "REBOUNDx will overwrite sim->additional_forces and sim->post_timestep_modifications.  If you are using custom functions you wrote and setting them through these function pointers, this won't work.  You should instead add your custom effects through REBOUNDx.  See http://reboundx.readthedocs.org/en/latest/c_examples.html#adding-custom-post-timestep-modifications-and-forces for a tutorial.");
    }
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
        free(current->param->contents);
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
    
    effect->hash = reb_hash(name);
    effect->ap = NULL;
    effect->force = NULL;
    effect->ptm = NULL;
    effect->rebx = rebx;
    
    effect->name = malloc(strlen(name) + 1); // +1 for \0 at end
    if (effect->name != NULL){
        strcpy(effect->name, name);
    }
    
    effect->next = rebx->effects;
    rebx->effects = effect;
   
    struct reb_particle* p = (struct reb_particle*)effect;
    p->ap = (void*)REBX_OBJECT_TYPE_EFFECT;
    return effect;
}

struct rebx_effect* rebx_add(struct rebx_extras* rebx, const char* name){
    struct rebx_effect* effect = rebx_add_effect(rebx, name);
    struct reb_simulation* sim = rebx->sim;

    if (effect->hash == reb_hash("modify_orbits_direct")){
        effect->ptm = rebx_modify_orbits_direct;
    }
    else if (effect->hash == reb_hash("modify_orbits_forces")){
        sim->force_is_velocity_dependent = 1;
        effect->force = rebx_modify_orbits_forces;
    }
    else if(effect->hash == reb_hash("gr")){
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
    else if (effect->hash == reb_hash("radiation_forces")){
        sim->force_is_velocity_dependent = 1;
        effect->force = rebx_radiation_forces;
    }
    else if (effect->hash == reb_hash("tides_precession")){
        effect->force = rebx_tides_precession;
    }
    else if (effect->hash == reb_hash("central_force")){
        effect->force = rebx_central_force;
    }
    else if (effect->hash == reb_hash("track_min_distance")){
        effect->ptm = rebx_track_min_distance;
    }
    else if (effect->hash == reb_hash("moon_quadrupole_laskar")){
        effect->force = rebx_moon_quadrupole_laskar;
    }
    else{
        char str[100]; 
        sprintf(str, "Effect '%s' passed to rebx_add not found.\n", name);
        reb_error(sim, str);
    }
    
    // Make sure function pointers in simulation are set, but only set them if there's at least one REBOUNDx effect using them.
    if(effect->force){
        sim->additional_forces = rebx_forces;
    }
    if(effect->ptm){
        sim->post_timestep_modifications = rebx_post_timestep_modifications;
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
    rebx->sim->additional_forces = rebx_forces;
    return effect;
}

struct rebx_effect* rebx_add_custom_post_timestep_modification(struct rebx_extras* rebx, const char* name, void (*custom_ptm)(struct reb_simulation* const sim, struct rebx_effect* const effect)){
    struct rebx_effect* effect = rebx_add_effect(rebx, name);
    effect->ptm = custom_ptm;
    rebx->sim->post_timestep_modifications = rebx_post_timestep_modifications;
    return effect;
}
    
void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param){
    struct rebx_param_to_be_freed* newparam = malloc(sizeof(*newparam));
    newparam->param = param;

    newparam->next = rebx->params_to_be_freed;
    rebx->params_to_be_freed = newparam;
}

/*********************************************************************************
 Internal functions for dealing with parameters
 ********************************************************************************/

static enum rebx_object_type rebx_get_object_type(const void* const object){
    struct reb_particle* p = (struct reb_particle*)object;
    if (p->ap == (void*)(intptr_t)REBX_OBJECT_TYPE_EFFECT){    // In add_effect we cast effect to a particle and set p->ap to REBX_OBJECT_TYPE_EFFECT
        return REBX_OBJECT_TYPE_EFFECT;
    }
    else{
        return REBX_OBJECT_TYPE_PARTICLE;
    }
}


// get simulation pointer from an object
static struct reb_simulation* rebx_get_sim(const void* const object){
    switch(rebx_get_object_type(object)){
        case REBX_OBJECT_TYPE_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            return effect->rebx->sim;
        }
        case REBX_OBJECT_TYPE_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            return p->sim;
        }
    }
    return NULL;
} 

struct rebx_param* rebx_create_param(){
    struct rebx_param* newparam = malloc(sizeof(*newparam));
    int collision=0;
    for(int j=REBX_OBJECT_TYPE_EFFECT; j<=REBX_OBJECT_TYPE_PARTICLE; j++){
        /* need this cast to avoid warnings (we are converting 32 bit int enum to 64 bit pointer). I think behavior
         * is implemetation defined, but this is OK because we are always checking the ap pointer in this way, so 
         * just need that whatever implementation defined result we get from this comparison is false for the address
         * of newparam, which will be the head node in the linked list at object->ap*/
        if(newparam == (void*)(intptr_t)j){ 
            collision=1;
        }
    }
    
    if (collision){
        free(newparam);
        struct rebx_param* address[5] = {NULL};
        int i;
        for(i=0; i<=5; i++){
            address[i] = malloc(sizeof(*newparam));
            collision = 0;
            for(int j=REBX_OBJECT_TYPE_EFFECT; j<=REBX_OBJECT_TYPE_PARTICLE; j++){
                if(address[i] == (void*)(intptr_t)j){
                    collision=1;
                }
            }
            if (collision == 0){
                newparam = address[i];
                break;
            }
        }

        if (i==5){
            fprintf(stderr, "REBOUNDx Error:  Can't allocate valid memory for parameter.\n");
            return NULL;
        }
        for(int j=0;j<i;j++){
            free(address[j]);
        }
    }
    return newparam;
}

/*********************************************************************************
 User interface for parameters
 ********************************************************************************/
int rebx_remove_param(const void* const object, const char* const param_name){
    // TODO free memory for deleted node
    uint32_t hash = reb_hash(param_name);
    struct rebx_param* current;
    switch(rebx_get_object_type(object)){
        case REBX_OBJECT_TYPE_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            current = effect->ap;
            if(current->hash == hash){
                effect->ap = current->next;
                return 1;
            }
            break;
        }
        case REBX_OBJECT_TYPE_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            current = p->ap;
            if(current->hash == hash){
                p->ap = current->next;
                return 1;
            }
            break;
        }
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

size_t rebx_sizeof(enum rebx_param_type param_type){
    switch(param_type){
        case REBX_TYPE_DOUBLE:
        {
            return sizeof(double);
        }
        case REBX_TYPE_INT:
        {
            return sizeof(int);
        }
        case REBX_TYPE_UINT32:
        {
            return sizeof(uint32_t);
        }
        case REBX_TYPE_ORBIT:
        {
            return sizeof(struct reb_orbit);
        }
        case REBX_TYPE_LONGLONG:
        {
            return sizeof(long long);
        }
    }
    return 0; // type not found
}

struct rebx_param* rebx_attach_param_node(void* const object, struct rebx_param* param){
    void* ptr = rebx_get_param(object, param->name);
    if (ptr != NULL){
        char str[300];
        sprintf(str, "REBOUNDx Error: Parameter '%s' passed to rebx_add_param already exists.\n", param->name);
        reb_error(rebx_get_sim(object), str);
        return NULL;
    }
    
    switch(rebx_get_object_type(object)){
        case REBX_OBJECT_TYPE_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            param->next = effect->ap;
            effect->ap = param;
            break;
        }
        case REBX_OBJECT_TYPE_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            param->next = p->ap;
            p->ap = param;
            break;
        }
    }
    
    return param;
}

struct rebx_param* rebx_add_param_node(void* const object, const char* const param_name, enum rebx_param_type param_type, const int ndim, const int* const shape){
    struct rebx_param* newparam = rebx_create_param();
    newparam->name = malloc(strlen(param_name) + 1); // +1 for \0 at end
    if (newparam->name != NULL){
        strcpy(newparam->name, param_name);
    }
    newparam->hash = reb_hash(param_name);
    newparam->param_type = param_type;
    newparam->python_type = -1; // not used by C
    newparam->ndim = ndim;
    newparam->shape = NULL;
    newparam->size = 1;
    if (ndim > 0){
	    size_t shapesize = sizeof(int)*ndim;
	    newparam->shape = malloc(shapesize);
        newparam->strides = malloc(shapesize);
	    memcpy(newparam->shape, shape, shapesize);
	    for(int i=ndim-1;i>=0;i--){ // going backward allows us to calculate strides at the same time
            newparam->strides[i] = newparam->size;  // stride[i] is equal to the product of the shapes for all indices > i
		    newparam->size *= shape[i];
	    }
    }
    size_t element_size = rebx_sizeof(param_type);
    if (element_size){
        newparam->contents = malloc(element_size*newparam->size); // newparam->size = number of elements in array (1 if scalar)
    }

    else{
        char str[300];
        sprintf(str, "REBOUNDx Error: Parameter type '%d' passed to rebx_add_param_node not supported.\n", param_type);
        reb_error(rebx_get_sim(object), str);
    }

    newparam = rebx_attach_param_node(object, newparam);
    return newparam;
}

void* rebx_add_param_array(void* const object, const char* const param_name, enum rebx_param_type param_type, const int ndim, const int* const shape){
    return rebx_add_param_node(object, param_name, param_type, ndim, shape)->contents;
}

void* rebx_add_param(void* const object, const char* const param_name, enum rebx_param_type param_type){
	int ndim=0;
	return rebx_add_param_array(object, param_name, param_type, ndim, NULL);
}

void* rebx_get_param(const void* const object, const char* const param_name){
    struct rebx_param* node = rebx_get_param_node(object, param_name);
    if (node == NULL){
        return NULL;
    }
    return node->contents;
}


struct rebx_param* rebx_get_param_node(const void* const object, const char* const param_name){
    struct rebx_param* current;
    switch(rebx_get_object_type(object)){
        case REBX_OBJECT_TYPE_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            current = effect->ap;
            break;
        }
        case REBX_OBJECT_TYPE_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            current = p->ap;
            break;
        }
    }
    
    uint32_t hash = reb_hash(param_name);
    while(current != NULL){
        if(current->hash == hash){
            return current; 
        }
        current = current->next;
    }
   
    if (current == NULL){   // param_name not found.  Return immediately.
        return NULL;
    }

    return current;
}

void* rebx_get_param_check(const void* const object, const char* const param_name, enum rebx_param_type param_type){
    struct rebx_param* node = rebx_get_param_node(object, param_name);
    if (node == NULL){
        return NULL;
    }
    
    if (node->param_type != param_type){
        char str[300];
        sprintf(str, "REBOUNDx Error: Parameter '%s' passed to rebx_get_param_check was found but was of wrong type.  See documentation for your particular effect.  In python, you might need to add a dot at the end of the number when assigning a parameter that REBOUNDx expects as a float.\n", param_name);
        reb_error(rebx_get_sim(object), str);
        return NULL;
    }

    return node->contents;
}

struct rebx_effect* rebx_get_effect(struct rebx_extras* const rebx, const char* const effect_name){
    struct rebx_effect* current = rebx->effects;
    uint32_t hash = reb_hash(effect_name);
    while(current != NULL){
        if(current->hash == hash){
            return current;
        }
        current = current->next;
    }
    
    if (current == NULL){   // effect_name not found.  Return immediately.
        return NULL;
    }
    
    return current;
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
