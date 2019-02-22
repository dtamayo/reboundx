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
#include "linkedlist.h"

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

void rebx_register_default_params(struct rebx_extras* rebx){
    rebx_register_param(rebx, "c", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "gr_source", REBX_TYPE_INT);
    rebx_register_param(rebx, "max_iterations", REBX_TYPE_INT);
    rebx_register_param(rebx, "tau_mass", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "index", REBX_TYPE_INT);
    rebx_register_param(rebx, "force", REBX_TYPE_FORCE);
    rebx_register_param(rebx, "particle", REBX_TYPE_POINTER);
}

void rebx_register_param(struct rebx_extras* const rebx, const char* name, enum rebx_param_type type){
    
    enum rebx_param_type reg_type = rebx_get_type(rebx, name);
    
    if (reg_type != REBX_TYPE_NONE){
        char str[300];
        sprintf(str, "REBOUNDx Error: Parameter name '%s' already in registered list. Cannot add duplicates.\n", name);
        reb_error(rebx->sim, str);
        return;
    }
    
    // Create new entry. These are just rebx_param structs without value populated
    struct rebx_param* param = rebx_create_param(rebx, name, type);
    if (param == NULL){
        return;
    }
    int success = rebx_add_param(rebx, &rebx->registered_params, param);
    if(!success){
        rebx_free_param(param);
    }
    
    return;
}

struct rebx_extras* rebx_attach(struct reb_simulation* sim){  // reboundx.h
    struct rebx_extras* rebx = malloc(sizeof(*rebx));
    rebx_initialize(sim, rebx);
    rebx_register_default_params(rebx);
    return rebx;
}

void rebx_detach(struct reb_simulation* sim){
    sim->additional_forces = NULL;
    sim->pre_timestep_modifications = NULL;
    sim->post_timestep_modifications = NULL;
    sim->free_particle_ap = NULL;
}

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx){
    sim->extras = rebx;
    rebx->sim = sim;
	//rebx->params_to_be_freed = NULL;
	rebx->additional_forces = NULL;
    rebx->pre_timestep_modifications=NULL;
    rebx->post_timestep_modifications=NULL;
    rebx->allocated_forces=NULL;
    rebx->allocated_operators=NULL;
    rebx->registered_params=NULL;
    rebx->integrator = REBX_INTEGRATOR_IMPLICIT_MIDPOINT;
    
    if(sim->additional_forces || sim->pre_timestep_modifications || sim->post_timestep_modifications){
        reb_warning(sim, "REBOUNDx overwrites sim->additional_forces, sim->pre_timestep_modifications and sim->post_timestep_modifications.  If you want to use REBOUNDx together with your own custom functions that use these callbacks, you should add them through REBOUNDx.  See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for a tutorial.");
    }
    // Have to set all the following at initialization since we can't know which will be needed ahead of time
    sim->additional_forces = rebx_additional_forces;
    sim->pre_timestep_modifications = rebx_pre_timestep_modifications;
    sim->post_timestep_modifications = rebx_post_timestep_modifications;
    sim->free_particle_ap = rebx_free_particle_ap;
}

void rebx_free(struct rebx_extras* rebx){
    rebx_free_pointers(rebx);
    free(rebx);
}

/**********************************************
 User Interface for adding forces and operators
 *********************************************/

struct rebx_force* rebx_create_force(struct rebx_extras* const rebx, const char* name){
    struct rebx_force* force = rebx_malloc(rebx, sizeof(*force));
    if (force == NULL){
        return NULL;
    }
    force->ap = NULL;
    force->sim = rebx->sim;
    force->force_type = REBX_FORCE_NONE;
    force->update_accelerations = NULL;
    force->name = NULL;
    if(name != NULL)
    {
        force->name = rebx_malloc(rebx, strlen(name) + 1); // +1 for \0 at end
        if (force->name == NULL){
            rebx_free_force(force);
            return NULL;
        }
        else{
            strcpy(force->name, name);
        }
    }
    
    // Add force to allocated_forces list for later freeing
    struct rebx_node* node = rebx_create_node(rebx);
    if (node == NULL){
        rebx_free_force(force);
        return NULL;
    }
    node->object = force;
    rebx_add_node(&rebx->allocated_forces, node);
    
    return force;
}

struct rebx_force* rebx_load_force(struct rebx_extras* const rebx, const char* name){
    struct rebx_force* force = NULL;
    
    if(strcmp(name, "gr") == 0){
        force = rebx_create_force(rebx, name);
        force->update_accelerations = rebx_gr;
        force->force_type = REBX_FORCE_VEL;
    }
    /*if (strcmp(name, "modify_orbits_forces") == 0){
     force = rebx_create_force(rebx, name);
     force->update_accelerations = rebx_modify_orbits_forces;
     force->force_type = REBX_FORCE_VEL;
     }
     
     else if (strcmp(name, "gr_full") == 0){
     force = rebx_create_force(rebx, name);
     force->update_accelerations = rebx_gr_full;
     force->force_type = REBX_FORCE_VEL;
     }
     else if (strcmp(name, "gravitational_harmonics") == 0){
     force = rebx_create_force(rebx, name);
     force->update_accelerations = rebx_gravitational_harmonics;
     force->force_type = REBX_FORCE_POS;
     }*/
    /*
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
     */
    else{
        char str[300];
        sprintf(str, "REBOUNDx error: Force '%s' not found in REBOUNDx library.\n", name);
        reb_error(rebx->sim, str);
    }
    
    return force; // NULL if not found
}

struct rebx_operator* rebx_create_operator(struct rebx_extras* const rebx, const char* name){
    struct rebx_operator* operator = rebx_malloc(rebx, sizeof(*operator));
    if (operator == NULL){
        return NULL;
    }
    operator->ap = NULL;
    operator->sim = rebx->sim;
    operator->operator_type = REBX_OPERATOR_NONE;
    operator->step = NULL;
    operator->name = NULL;
    if(name != NULL){
        operator->name = rebx_malloc(rebx, strlen(name) + 1); // +1 for \0 at end
        if (operator->name == NULL){
            rebx_free_operator(operator);
            return NULL;
        }
        else{
            strcpy(operator->name, name);
        }
    }
    
    // Add operator to allocated_operators list for later freeing
    struct rebx_node* node = rebx_create_node(rebx);
    if (node == NULL){
        rebx_free_operator(operator);
        return 0;
    }
    node->object = operator;
    rebx_add_node(&rebx->allocated_operators, node);
    
    return operator;
}

struct rebx_operator* rebx_load_operator(struct rebx_extras* const rebx, const char* name){
    struct rebx_operator* operator = NULL;
    
    if (strcmp(name, "modify_mass") == 0){
        operator = rebx_create_operator(rebx, name);
        operator->step = rebx_modify_mass;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }
    /*else if (strcmp(name, "kepler") == 0){
        operator = rebx_create_operator(rebx, name);
        operator->step = rebx_kepler_step;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }
    else if (strcmp(name, "jump") == 0){
        operator = rebx_create_operator(rebx, name);
        operator->step = rebx_jump_step;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }
    else if (strcmp(name, "interaction") == 0){
        operator = rebx_create_operator(rebx, name);
        operator->step = rebx_interaction_step;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }
    else if (strcmp(name, "ias15") == 0){
        operator = rebx_create_operator(rebx, name);
        operator->step = rebx_ias15_step;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }
    
     else if (strcmp(name, "modify_orbits_direct") == 0){
     operator = rebx_create_operator(rebx, name);
     operator->step = rebx_modify_orbits_direct;
     operator->operator_type = REBX_OPERATOR_UPDATER;
     }
     */

    /*else if (hash == reb_hash("track_min_distance")){
     operator->step = rebx_track_min_distance;
     }*/
    else{
        char str[300];
        sprintf(str, "REBOUNDx error: Operator '%s' not found in REBOUNDx library.\n", name);
        reb_error(rebx->sim, str);
    }
    return operator; // NULL if not found
}

int rebx_add_force(struct rebx_extras* rebx, struct rebx_force* force){
    if (force == NULL){
        return 0;
    }
    
    if (force->update_accelerations == NULL){
        reb_error(rebx->sim, "REBOUNDx error: Need to set update_accelerations function pointer on force before calling rebx_add_force. See custom effects example.\n");
        return 0;
    }
    
    if (force->force_type == REBX_FORCE_NONE){
        reb_error(rebx->sim, "REBOUNDx error: Need to set force_type field on force before calling rebx_add_force. See custom effects example.\n");
        return 0;
    }
    
    if (force->force_type == REBX_FORCE_VEL){
        rebx->sim->force_is_velocity_dependent = 1;
    }
    
    // Could add logic based on different integrators
    struct rebx_node* node = rebx_create_node(rebx);
    if (node == NULL){
        return 0;
    }
    node->object = force;
    rebx_add_node(&rebx->additional_forces, node);
    
    return 1;
}

int rebx_add_operator_step(struct rebx_extras* rebx, struct rebx_operator* operator, const double dt_fraction, enum rebx_timing timing){
    if (operator->step == NULL){
        reb_error(rebx->sim, "REBOUNDx error: Need to set step function pointer on operator before adding to simulation. See custom effects example.\n");
        return 0;
    }
    
    if (operator->operator_type == REBX_OPERATOR_NONE){
        reb_error(rebx->sim, "REBOUNDx error: Need to set operator_type field on operator before adding to simulation. See custom effects example.\n");
        return 0;
    }
    
    struct rebx_step* step = rebx_malloc(rebx, sizeof(*step));
    if(step == NULL){
        return 0;
    }
    step->operator = operator;
    step->dt_fraction = dt_fraction;
    
    struct rebx_node* node = rebx_create_node(rebx);
    if (node == NULL){
        return 0;
    }
    node->object = step;
    if (timing == REBX_TIMING_PRE){
        rebx_add_node(&rebx->pre_timestep_modifications, node);
        return 1;
    }
    if (timing == REBX_TIMING_POST){
        rebx_add_node(&rebx->post_timestep_modifications, node);
        return 1;
    }
    return 0;
}
    
int rebx_add_operator(struct rebx_extras* rebx, struct rebx_operator* operator){
    if (operator == NULL){
        return 0;
    }
    
    struct reb_simulation* const sim = rebx->sim;
    double dt_fraction;
    if (operator->operator_type == REBX_OPERATOR_RECORDER){
        // Doesn't alter state. Add once after timestep.
        dt_fraction = 1.;
        int success = rebx_add_operator_step(rebx, operator, dt_fraction, REBX_TIMING_POST);
        return success;
    }
    
    switch(sim->integrator){
        case REB_INTEGRATOR_IAS15:
        // don't add pre-timestep b/c don't know what IAS will choose as dt
        {
            dt_fraction = 1.;
            int success = rebx_add_operator_step(rebx, operator, dt_fraction, REBX_TIMING_POST);
            return success;
        }
        case REB_INTEGRATOR_WHFAST: // half step pre and post
        {
            dt_fraction = 1./2.;
            int success1 = rebx_add_operator_step(rebx, operator, dt_fraction, REBX_TIMING_PRE);
            int success2 = rebx_add_operator_step(rebx, operator, dt_fraction, REBX_TIMING_POST);
            return (success1 && success2);
        }
        case REB_INTEGRATOR_MERCURIUS: // half step pre and post
        {
            if (operator->operator_type == REBX_OPERATOR_UPDATER){
                reb_error(sim, "REBOUNDx Error: Operators that change particle states are not supported with Mercurius.\n");
                return 0;
            }
        }
    }
    return 0; // didn't reach a successful outcome
}

/*******************************************************************
 User interface for setting parameter values
 *******************************************************************/

int rebx_set_param_pointer(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, void* val){
    if (apptr == NULL){
        reb_error(rebx->sim, "REBOUNDx Error: Passed NULL apptr to rebx_add_param. See examples.\n");
        return 0;
    }
    
    enum rebx_param_type type = rebx_get_type(rebx, param_name);
    if (type == REBX_TYPE_NONE){
        char str[300];
        sprintf(str, "REBOUNDx Error: Need to register parameter name '%s' before using it. See examples.\n", param_name);
        reb_error(rebx->sim, str);
        return 0;
    }
    
    // Check whether it already exists in linked list
    struct rebx_param* param = rebx_get_param_struct(rebx, *apptr, param_name);
    
    if(param == NULL){
        param = rebx_create_param(rebx, param_name, type);
        if (param == NULL){ // adding new param failed
            return 0;
        }
        int success = rebx_add_param(rebx, apptr, param);
        if(!success){
            rebx_free_param(param);
            return 0;
        }
    }
    
    param->value = val;
    return 1;
}

int rebx_set_param_double(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, double val){
    if (apptr == NULL){
        reb_error(rebx->sim, "REBOUNDx Error: Passed NULL apptr to rebx_set_param_double. See examples.\n");
        return 0;
    }
    
    enum rebx_param_type type = rebx_get_type(rebx, param_name);
    if (type == REBX_TYPE_NONE){
        char str[300];
        sprintf(str, "REBOUNDx Error: Need to register parameter name '%s' before using it. See examples.\n", param_name);
        reb_error(rebx->sim, str);
        return 0;
    }
    
    // Check whether param already exists in linked list
    struct rebx_param* param = rebx_get_param_struct(rebx, *apptr, param_name);
    
    if(param == NULL){      // Make new param and allocate memory
        param = rebx_create_param(rebx, param_name, type);
        if (param == NULL){ // adding new param failed
            return 0;
        }
        int success = rebx_add_param(rebx, apptr, param);
        if(!success){
            rebx_free_param(param);
            return 0;
        }
        param->value = rebx_malloc(rebx, sizeof(double));
    }
    // Update new or existing param value
    double* valptr = param->value;
    *valptr = val;
    
    return 1;
}

int rebx_set_param_int(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, int val){
    if (apptr == NULL){
        reb_error(rebx->sim, "REBOUNDx Error: Passed NULL apptr to rebx_add_param. See examples.\n");
        return 0;
    }
    
    enum rebx_param_type type = rebx_get_type(rebx, param_name);
    if (type == REBX_TYPE_NONE){
        char str[300];
        sprintf(str, "REBOUNDx Error: Need to register parameter name '%s' before using it. See examples.\n", param_name);
        reb_error(rebx->sim, str);
        return 0;
    }
    
    // Check whether it already exists in linked list
    struct rebx_param* param = rebx_get_param_struct(rebx, *apptr, param_name);
    
    if(param == NULL){      // Make new param and allocate memory
        param = rebx_create_param(rebx, param_name, type);
        if (param == NULL){ // adding new param failed
            return 0;
        }
        int success = rebx_add_param(rebx, apptr, param);
        if(!success){
            rebx_free_param(param);
            return 0;
        }
        param->value = rebx_malloc(rebx, sizeof(int));
    }
    int* valptr = param->value;
    *valptr = val;
    
    return 1;
}

/*******************************************************************
 User interface for getting REBOUNDx objects and parameters
 *******************************************************************/

struct rebx_param* rebx_get_param_struct(struct rebx_extras* rebx, struct rebx_node* ap, const char* const param_name){
    struct rebx_node* current = ap;
    while(current != NULL){
        struct rebx_param* param = current->object;
        if(strcmp(param->name, param_name) == 0){
            return param;
        }
        current = current->next;
    }
    
    return NULL;   // name not found. Don't want warnings for optional parameters so don't reb_error
}

void* rebx_get_param(struct rebx_extras* rebx, struct rebx_node* ap, const char* const param_name){
    struct rebx_param* param = rebx_get_param_struct(rebx, ap, param_name);
    if (param == NULL){
        return NULL;
    }
    else{
        return param->value;
    }
}

struct rebx_force* rebx_get_force(struct rebx_extras* const rebx, const char* const name){
    struct rebx_node* current = rebx->allocated_forces;
    while(current != NULL){
        struct rebx_force* force = current->object;
        if(strcmp(force->name, name) == 0){
            return force;
        }
        current = current->next;
    }
    
    return NULL;
}

struct rebx_operator* rebx_get_operator(struct rebx_extras* const rebx, const char* const name){
    struct rebx_node* current = rebx->allocated_operators;
    while(current != NULL){
        struct rebx_operator* operator = current->object;
        if(strcmp(operator->name, name) == 0){
            return operator;
        }
        current = current->next;
    }
    
    return NULL;
}

/*******************************************************************
 User interface for removing REBOUNDx objects
 *******************************************************************/

int rebx_remove_force(struct rebx_extras* rebx, struct rebx_force* force){
    int allocated = rebx_remove_node(&rebx->allocated_forces, force);
    if(allocated){
        rebx_free_force(force);
    }
    // success only cares about removal from add_forces that affects sim
    int success = rebx_remove_node(&rebx->additional_forces, force);
    return success;
}

// Remove all steps in head pointer that have the passed operator in them.
// Success = 1 if at least one removed. Need separate logic since operator
// is nested inside step
static int rebx_remove_step_node(struct rebx_node** head, struct rebx_operator* operator){
    if (*head == NULL){
        return 0;
    }
    
    struct rebx_node* current = *head;
    struct rebx_step* step = current->object;
    if(step->operator == operator){ // edge case where step is first in list
        *head = current->next;
        rebx_free_step(step);
        free(current);
        return 1;
    }
    
    struct rebx_node* prev = current;
    current = current->next;
    while (current != NULL){
        step = current->object;
        if(step->operator == operator){
            prev->next = current->next;
            rebx_free_step(step);
            free(current);
            return 1;
        }
        prev = current;
        current = current->next;
    }
    return 0;
}

int rebx_remove_operator(struct rebx_extras* rebx, struct rebx_operator* operator){
    int allocated = rebx_remove_node(&rebx->allocated_operators, operator);
    if(allocated){
        rebx_free_operator(operator);
        
    }
    
    // success only cares about removal from lists that actually do
    // something to sim below. Success if EITHER one successful.
    int success = 0;
    int keep_searching = 1;
    while(keep_searching){ // keep searching while steps are found
        keep_searching = rebx_remove_step_node(&rebx->pre_timestep_modifications, operator);
        if (keep_searching == 1){ // success if at least one step found
            success = 1;
        }
    }

    keep_searching = 1;
    while(keep_searching){ // keep searching while steps are found
        keep_searching = rebx_remove_step_node(&rebx->post_timestep_modifications, operator);
        if (keep_searching == 1){ // success if at least one step found
            success = 1;
        }
    }
    
    return success;
}

/***************************************************************
 * Internal Memory Handling Routines
 ******************************************************************/

void* rebx_malloc(struct rebx_extras* const rebx, size_t memsize){
    void* ptr = malloc(memsize);
    if (ptr == NULL && memsize>0){
        reb_error(rebx->sim, "REBOUNDx Error: Could not allocate memory.\n");
        return NULL;
    }
    
    return ptr;
}

void rebx_free_param(struct rebx_param* param){
    if(param->name){
        free(param->name);
    }
    // Don't free pointers to structs
    if(param->type == REBX_TYPE_INT || param->type == REBX_TYPE_DOUBLE){
        if(param->value){
            free(param->value);
        }
    }
    free(param);
}

void rebx_free_ap(struct rebx_node** ap){
    struct rebx_node* current = *ap;
    struct rebx_node* next;
    while (current != NULL){
        next = current->next;
        rebx_free_param(current->object);
        free(current);
        current = next;
    }
}

void rebx_free_particle_ap(struct reb_particle* p){
    rebx_free_ap(&p->ap);
}

void rebx_free_force(struct rebx_force* force){
    if(force->name){
        free(force->name);
    }
    rebx_free_ap(&force->ap);
    free(force);
}

void rebx_free_operator(struct rebx_operator* operator){
    if(operator->name){
        free(operator->name);
    }
    rebx_free_ap(&operator->ap);
    free(operator);
}

void rebx_free_step(struct rebx_step* step){
    free(step);
}

void rebx_free_reg_param(struct rebx_param* param){
    if(param->name){
        free(param->name);
    }
    free(param);
}

void rebx_free_pointers(struct rebx_extras* rebx){
    struct rebx_node* current;
    struct rebx_node* next;
    
    current = rebx->allocated_forces;
    while (current != NULL){
        next = current->next;
        rebx_free_force(current->object);
        free(current);
        current = next;
    }
    
    current = rebx->allocated_operators;
    while (current != NULL){
        next = current->next;
        rebx_free_operator(current->object);
        free(current);
        current = next;
    }
    
    current = rebx->additional_forces;
    while (current != NULL){
        next = current->next;
        free(current);
        current = next;
    }
    
    
    current = rebx->pre_timestep_modifications;
    while (current != NULL){
        next = current->next;
        rebx_free_step(current->object);
        free(current);
        current = next;
    }
    
    current = rebx->post_timestep_modifications;
    while (current != NULL){
        next = current->next;
        rebx_free_step(current->object);
        free(current);
        current = next;
    }
    
    current = rebx->registered_params;
    while (current != NULL){
        next = current->next;
        rebx_free_reg_param(current->object);
        free(current);
        current = next;
    }
}

/**********************************************
 Internal Functions executing forces & ptm each timestep
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
    struct rebx_node* current = rebx->additional_forces;
    while(current != NULL){
        /*if(sim->force_is_velocity_dependent && sim->integrator==REB_INTEGRATOR_WHFAST){
         reb_warning(sim, "REBOUNDx: Passing a velocity-dependent force to WHFAST. Need to apply as an operator.");
         }*/
        struct rebx_force* force = current->object;
        const double N = sim->N - sim->N_var;
        force->update_accelerations(sim, force, sim->particles, N);
        current = current->next;
    }
}

void rebx_pre_timestep_modifications(struct reb_simulation* sim){
    struct rebx_extras* rebx = sim->extras;
    struct rebx_node* current = rebx->pre_timestep_modifications;
    const double dt = sim->dt;
    
    while(current != NULL){
        if(sim->integrator==REB_INTEGRATOR_IAS15 && sim->ri_ias15.epsilon != 0){
            reb_warning(sim, "REBOUNDx: Can't use pre-timestep modifications with adaptive timesteps (IAS15).");
        }
        struct rebx_step* step = current->object;
        struct rebx_operator* operator = step->operator;
        operator->step(sim, operator, dt*step->dt_fraction);
        current = current->next;
    }
}

void rebx_post_timestep_modifications(struct reb_simulation* sim){
    struct rebx_extras* rebx = sim->extras;
    struct rebx_node* current = rebx->post_timestep_modifications;
    const double dt = sim->dt_last_done;
    
    while(current != NULL){
        struct rebx_step* step = current->object;
        struct rebx_operator* operator = step->operator;
        operator->step(sim, operator, dt*step->dt_fraction);
        current = current->next;
    }
}

/****************************************************************
 Internal functions for dealing with parameters
 ****************************************************************/

struct rebx_node* rebx_create_node(struct rebx_extras* rebx){
    struct rebx_node* node = rebx_malloc(rebx, sizeof(*node));
    if (node == NULL){
        return NULL;
    }
    node->object = NULL;
    node->next = NULL;
    return node;
}

struct rebx_param* rebx_create_param(struct rebx_extras* rebx, const char* name, enum rebx_param_type type){
    // Allocate and initialize new param struct
    struct rebx_param* param = rebx_malloc(rebx, sizeof(*param));
    if (param == NULL){
        return NULL;
    }
    param->type = type;
    param->value = NULL;
    param->name = rebx_malloc(rebx, strlen(name) + 1); // +1 for \0 at end
    if (param->name == NULL){
        return NULL;
    }
    else{
        strcpy(param->name, name);
    }
    
    return param;
}

int rebx_add_param(struct rebx_extras* const rebx, struct rebx_node** apptr, struct rebx_param* param){
    struct rebx_node* node = rebx_create_node(rebx);
    if (node == NULL){
        return 0;
    }
    node->object = param;
    rebx_add_node(apptr, node);
    return 1;
}

// needed from Python
enum rebx_param_type rebx_get_type(struct rebx_extras* rebx, const char* name){
    struct rebx_param* param = rebx_get_param_struct(rebx, rebx->registered_params, name);
    
    if (param == NULL){ // param not found
        return REBX_TYPE_NONE;
    }
    
    return param->type;
}

size_t rebx_sizeof(struct rebx_extras* rebx, enum rebx_param_type type){
    switch(type){
        case REBX_TYPE_DOUBLE:
        {
            return sizeof(double);
        }
        case REBX_TYPE_INT:
        {
            return sizeof(int);
        }
        case REBX_TYPE_FORCE:
        {
            return sizeof(struct rebx_force);
        }
        case REBX_TYPE_POINTER:
        {
            return 0;
        }
        case REBX_TYPE_NONE:
        {
            reb_error(rebx->sim, "REBOUNDx Error: Parameter name passed to rebx_sizeof was not registered. This should not happen. Please open issue on github.com/dtamayo/reboundx.\n");
            return 0;
        }
        default:
        {
            reb_error(rebx->sim, "REBOUNDx Error: Need to add new param type to switch statement in rebx_sizeof. Please open issue on github.com/dtamayo/reboundx.\n");
            return 0;
        }
    }
}
