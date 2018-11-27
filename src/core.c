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
#include <float.h>
#include "core.h"
#include "rebound.h"

#define STRINGIFY(s) str(s)
#define str(s) #s

const char* rebx_build_str = __DATE__ " " __TIME__; // Date and time build string. 
const char* rebx_version_str = "2.19.3";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* rebx_githash_str = STRINGIFY(REBXGITHASH);             // This line gets updated automatically. Do not edit manually.

void rebx_integrate(struct reb_simulation* const sim, const double dt, struct rebx_effect* const effect){
    /*if (effect->force == NULL){
        char str[300];
        sprintf(str, "REBOUNDx Error: rebx_integrate called with non-force effect '%s'.\n", effect->name);
        reb_error(sim, str);
    }
    struct rebx_extras* rebx = sim->extras;
    rebx_reset_accelerations(sim->particles, sim->N);
    
    switch(rebx->integrator){
        case REBX_INTEGRATOR_IMPLICIT_MIDPOINT:
            rebx_integrator_implicit_midpoint_integrate(sim, dt, effect);
            break;
        case REBX_INTEGRATOR_RK2:
            rebx_integrator_rk2_integrate(sim, dt, effect);
            break;
        case REBX_INTEGRATOR_RK4:
            rebx_integrator_rk4_integrate(sim, dt, effect);
            break;
        case REBX_INTEGRATOR_EULER:
            rebx_integrator_euler_integrate(sim, dt, effect);
            break;
        case REBX_INTEGRATOR_NONE:
            break;
        default:
            break;
    }*/
}

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
	rebx->forces = NULL;
    rebx->pre_timestep_operators=NULL;
    rebx->post_timestep_operators=NULL;
    rebx->integrator = REBX_INTEGRATOR_IMPLICIT_MIDPOINT;
    
    if(sim->additional_forces || sim->pre_timestep_modifications || sim->post_timestep_modifications){
        reb_warning(sim, "REBOUNDx overwrites sim->additional_forces, sim->pre_timestep_modifications and sim->post_timestep_modifications.  If you want to use REBOUNDx together with your own custom functions that use these callbacks, you should add them through REBOUNDx.  See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for a tutorial.");
    }
    // Have to set all the following at initialization since we can't know
    // which will be needed from added effects. User could set force_as_operator after the fact.
    sim->additional_forces = rebx_additional_forces;
    sim->pre_timestep_modifications = rebx_pre_timestep_modifications;
    sim->post_timestep_modifications = rebx_post_timestep_modifications;
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
    /*struct rebx_effect* current = rebx->effects;
    struct rebx_effect* temp_next;

    while(current != NULL){
        temp_next = current->next;
        free(current);
        current = temp_next;
    }*/
}

/**********************************************
 Functions executing forces & ptm each timestep
 *********************************************/

void rebx_reset_accelerations(struct reb_particle* const ps, const int N){
    for(int i=0; i<N; i++){
        ps[i].ax = 0.;
        ps[i].ay = 0.;
        ps[i].az = 0.;
    }
}

void rebx_additional_forces(struct reb_simulation* sim){
	struct rebx_extras* rebx = sim->extras;
	struct rebx_force* current = rebx->forces;
	while(current != NULL){
        /*if(sim->force_is_velocity_dependent && sim->integrator==REB_INTEGRATOR_WHFAST){
            reb_warning(sim, "REBOUNDx: Passing a velocity-dependent force to WHFAST. Need to apply as an operator.");
        }*/
        const double N = sim->N - sim->N_var;
		current->update_accelerations(sim, current->effect, sim->particles, N);
		current = current->next;
	}
}

void rebx_pre_timestep_modifications(struct reb_simulation* sim){
    struct rebx_extras* rebx = sim->extras;
    struct rebx_operator* current = rebx->pre_timestep_operators;
    
    while(current != NULL){
        if(sim->integrator==REB_INTEGRATOR_IAS15 && sim->ri_ias15.epsilon != 0){
            reb_warning(sim, "REBOUNDx: Can't use second order scheme with adaptive timesteps (IAS15). Must use operator_order = 1 or apply as force to get sensible results.");
        }
        current->step(sim, current->effect, current->dt);
        current = current->next;
    }
}

void rebx_post_timestep_modifications(struct reb_simulation* sim){
	struct rebx_extras* rebx = sim->extras;
	struct rebx_operator* current = rebx->post_timestep_operators;
   
    //const double dt = sim->dt_last_done;
    // first do the 2nd order operators for half a timestep, in reverse order
    
    while(current != NULL){
        current->step(sim, current->effect, current->dt);
        current = current->prev;
    }
}

/**********************************************
 Adders for linked lists in extras structure
 *********************************************/

static enum rebx_effect_type rebx_get_effect_type(struct rebx_extras* rebx, const char* name){
    uint32_t hash = reb_hash(name);
    if (hash == reb_hash("modify_orbits_forces")){
        return REBX_EFFECT_FORCE_VEL;
    }
    else if(hash == reb_hash("gr")){
        return REBX_EFFECT_FORCE_VEL;
    }
    else if (hash == reb_hash("modify_orbits_direct")){
        return REBX_EFFECT_OPERATOR_UPDATER;
    }
    else if (hash == reb_hash("modify_mass")){
        return REBX_EFFECT_OPERATOR_UPDATER;
    }
    else if (hash == reb_hash("track_min_distance")){
        return REBX_EFFECT_OPERATOR_RECORDER;
    }
    else{
        char str[100];
        sprintf(str, "REBOUNDx error: Effect '%s' not found in rebx_get_effect_type.\n", name);
        reb_error(rebx->sim, str);
        return REBX_EFFECT_NONE;
    }
}

static int rebx_set_force_fp(struct rebx_extras* rebx, struct rebx_force* force){
    uint32_t hash = reb_hash(force->effect->name);
    if (hash == reb_hash("modify_orbits_forces")){
        force->update_accelerations = rebx_modify_orbits_forces;
    }
    else if(hash == reb_hash("gr")){
        force->update_accelerations = rebx_gr;
    }
    /*else if (hash == reb_hash("gr_full")){
        sim->force_is_velocity_dependent = 1;
        update_accelerations = rebx_gr_full;
    }
    else if (hash == reb_hash("gr_potential")){
        update_accelerations = rebx_gr_potential;
    }
    else if (hash == reb_hash("radiation_forces")){
        sim->force_is_velocity_dependent = 1;
        update_accelerations = rebx_radiation_forces;
    }
    else if (hash == reb_hash("tides_precession")){
        update_accelerations = rebx_tides_precession;
    }
    else if (hash == reb_hash("central_force")){
        update_accelerations = rebx_central_force;
    }
    else if (hash == reb_hash("tides_synchronous_ecc_damping")){
        sim->force_is_velocity_dependent = 1;
        update_accelerations = rebx_tides_synchronous_ecc_damping;
    }
    else if (hash == reb_hash("gravitational_harmonics")){
        update_accelerations = rebx_gravitational_harmonics;
    }*/
    else{
        char str[100];
        sprintf(str, "REBOUNDx error: Effect '%s' not found in rebx_set_force_fp.\n", force->effect->name);
        reb_error(rebx->sim, str);
        return 0;
    }
    return 1;
}

static int rebx_set_operator_fp(struct rebx_extras* rebx, struct rebx_operator* operator){
    uint32_t hash = reb_hash(operator->effect->name);
    if (hash == reb_hash("modify_orbits_direct")){
        operator->step = rebx_modify_orbits_direct;
    }
    else if (hash == reb_hash("modify_mass")){
        operator->step = rebx_modify_mass;
    }
    else if (hash == reb_hash("track_min_distance")){
        operator->step = rebx_track_min_distance;
    }
    else{
        char str[100];
        sprintf(str, "REBOUNDx error: Effect '%s' not found in rebx_set_operator_fp.\n", operator->effect->name);
        reb_error(rebx->sim, str);
    }
    return 1;
}

struct rebx_effect* rebx_add_effect(struct rebx_extras* rebx, const char* name){
    struct rebx_effect* effect = malloc(sizeof(*effect));
    
    effect->hash = reb_hash(name);
    effect->ap = NULL;
    effect->rebx = rebx;
    
    effect->name = malloc(strlen(name) + 1); // +1 for \0 at end
    if (effect->name != NULL){
        strcpy(effect->name, name);
    }
    
    struct reb_particle* p = (struct reb_particle*)effect;
    p->ap = (void*)REBX_OBJECT_TYPE_EFFECT;
    return effect;
}

// put incompatible integrator warnings in add_force and add_operator
void rebx_add_force(struct rebx_extras* rebx, struct rebx_effect* effect){
    struct reb_simulation* sim = rebx->sim;
    if(effect == NULL){
        reb_error(sim, "REBOUNDx error: Effect passed to manual add function was NULL. Need to call rebx_add_effect first, or use rebx_add.\n");
        return;
    }
    struct rebx_force* force = malloc(sizeof(*force));
    force->effect = effect;
    
    if(rebx_set_force_fp(rebx, force)){
        if(rebx_get_effect_type(rebx, force->effect->name) == REBX_EFFECT_FORCE_VEL){
            sim->force_is_velocity_dependent = 1;
        }
        force->next = rebx->forces;
        force->prev = NULL;
        rebx->forces = force;
        if (force->next){
            force->next->prev = force;
        }
    }
    else{
        free(force);
    }
}

static struct rebx_operator* rebx_add_operator(struct rebx_extras* rebx, struct rebx_effect* effect, const double dt){
    struct reb_simulation* sim = rebx->sim;
    if(effect == NULL){
        reb_error(sim, "REBOUNDx error: Effect passed to manual add function was NULL. Need to call rebx_add_effect first, or use rebx_add.\n");
        return NULL;
    }
    struct rebx_operator* operator = malloc(sizeof(*operator));
    operator->dt = dt;
    operator->effect = effect;
    
    if(rebx_set_operator_fp(rebx, operator)){
        return operator;
    }
    else{
        free(operator);
        return NULL;
    }
}
    
void rebx_add_pre_timestep_operator(struct rebx_extras* rebx, struct rebx_effect* effect, const double dt){
    struct rebx_operator* operator = rebx_add_operator(rebx, effect, dt);
    if(operator != NULL){
        operator->next = rebx->pre_timestep_operators;
        operator->prev = NULL;
        rebx->pre_timestep_operators = operator;
        if (operator->next){
            operator->next->prev = operator;
        }
    }
}

void rebx_add_post_timestep_operator(struct rebx_extras* rebx, struct rebx_effect* effect, const double dt){
    struct rebx_operator* operator = rebx_add_operator(rebx, effect, dt);
    if(operator != NULL){
        operator->next = rebx->post_timestep_operators;
        operator->prev = NULL;
        rebx->post_timestep_operators = operator;
        if (operator->next){
            operator->next->prev = operator;
        }
    }
}

struct rebx_effect* rebx_add(struct rebx_extras* rebx, const char* name){
    struct rebx_effect* effect = rebx_add_effect(rebx, name);
    struct reb_simulation* sim = rebx->sim;
    enum rebx_effect_type type = rebx_get_effect_type(rebx, reb_hash(name));
    /*
    // put recorders to float up. Separate linked list for recorders?
    if(type == REBX_EFFECT_OPERATOR_RECORDER){
        rebx_add_post_timestep_operator(rebx, name, sim->dt);
    }
    
    if(sim->integrator == REB_INTEGRATOR_IAS15){
        if(type == REBX_EFFECT_FORCE_POS || type == REBX_EFFECT_FORCE_VEL){
            rebx_add_force(rebx, name);
        }
        else if(type == REBX_EFFECT_OPERATOR_UPDATER){
            if(sim->ri_ias15.epsilon != 0.){
                char str[100];
                sprintf(str, "REBOUNDx error: Operator '%s' not supported with IAS15 using a variable stepsize. Set sim.ri_ias15.epsilon=0 \n", name);
                reb_error(sim, str);
                return NULL;
            }
            else{
                rebx_add_pre_timestep_operator(rebx, name, sim->dt/2.);
                rebx_add_post_timestep_operator(rebx, name, sim->dt/2.);
            }
        }
    }
    else if(sim->integrator == REB_INTEGRATOR_WHFAST){
        if(type == REBX_EFFECT_FORCE_POS){
            rebx_add_force(rebx, name);
        }
        else if(type == REBX_EFFECT_FORCE_VEL){ //update to add rebx_integrator operator?
            rebx_add_force(rebx, name);
        }
        else if(type == REBX_EFFECT_OPERATOR_UPDATER){
            rebx_add_pre_timestep_operator(rebx, name, sim->dt/2.);
            rebx_add_post_timestep_operator(rebx, name, sim->dt/2.);
        }
    }
    else if(sim->integrator == REB_INTEGRATOR_MERCURIUS){
        if(type == REBX_EFFECT_FORCE_POS || type == REBX_EFFECT_FORCE_VEL){
            rebx_add_force(rebx, name);
        }
        else if(type == REBX_EFFECT_OPERATOR_UPDATER){
            char str[100];
            sprintf(str, "REBOUNDx error: Operator '%s' not supported with Mercurius integrator.\n", name);
            reb_error(sim, str);
            return NULL;
        }
    }
    */
    return NULL;
}

struct rebx_effect* rebx_add_custom_force(struct rebx_extras* rebx, const char* name, void (*custom_force)(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N), const int force_is_velocity_dependent){
    struct rebx_effect* effect = rebx_add_effect(rebx, name);
    /*effect->force = custom_force;
    struct reb_simulation* sim = rebx->sim;
    if(force_is_velocity_dependent){
        sim->force_is_velocity_dependent = 1;
    }*/
    return effect;
}

struct rebx_effect* rebx_add_custom_operator(struct rebx_extras* rebx, const char* name, void (*custom_operator)(struct reb_simulation* const sim, struct rebx_effect* const effect, const double dt, enum rebx_timing timing)){
    struct rebx_effect* effect = rebx_add_effect(rebx, name);
    //effect->operator = custom_operator;
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
    struct rebx_effect* current = rebx->forces->effect; // stopgap fix!
    /*uint32_t hash = reb_hash(effect_name);
    while(current != NULL){
        if(current->hash == hash){
            return current;
        }
        current = current->next;
    }
    
    if (current == NULL){   // effect_name not found.  Return immediately.
        return NULL;
    }
    */
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
