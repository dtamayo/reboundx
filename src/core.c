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
#define PRINT_MSG_2SX(ARG0, ARG1) fprintf(stderr, "%s : 0x%016llX\n", ARG0, (unsigned long long)ARG1)
void print_ptraddress(void** ptraddress){
    PRINT_MSG_2SX("Address of Pointer:", ptraddress);

}
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

static void rebx_register_params(struct rebx_extras* rebx){
    rebx_register_param(rebx, "c", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "max_iterations", REBX_TYPE_INT);
    rebx_register_param(rebx, "index", REBX_TYPE_INT);
    rebx_register_param(rebx, "force", REBX_TYPE_FORCE);
    rebx_register_param(rebx, "particle", REBX_TYPE_POINTER);
}

struct rebx_extras* rebx_init(struct reb_simulation* sim){  // reboundx.h
    struct rebx_extras* rebx = malloc(sizeof(*rebx));
    rebx_initialize(sim, rebx);
    return rebx;
}

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx){
    sim->extras = rebx;
    rebx->sim = sim;
	//rebx->params_to_be_freed = NULL;
	rebx->forces = NULL;
    rebx->pre_timestep_operators=NULL;
    rebx->post_timestep_operators=NULL;
    rebx->allocated_pointers=NULL;
    rebx->registered_params=NULL;
    rebx->integrator = REBX_INTEGRATOR_IMPLICIT_MIDPOINT;
    
    if(sim->additional_forces || sim->pre_timestep_modifications || sim->post_timestep_modifications){
        reb_warning(sim, "REBOUNDx overwrites sim->additional_forces, sim->pre_timestep_modifications and sim->post_timestep_modifications.  If you want to use REBOUNDx together with your own custom functions that use these callbacks, you should add them through REBOUNDx.  See https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Custom_Effects.ipynb for a tutorial.");
    }
    // Have to set all the following at initialization since we can't know
    // which will be needed from added effects. User could set force_as_operator after the fact.
    sim->additional_forces = rebx_additional_forces;
    sim->pre_timestep_modifications = rebx_pre_timestep_modifications;
    sim->post_timestep_modifications = rebx_post_timestep_modifications;
    
    rebx_register_params(rebx);
}

/*****************************
 Garbage Collection Routines
 ****************************/

void rebx_remove_from_simulation(struct reb_simulation* sim){
    sim->additional_forces = NULL;
    sim->post_timestep_modifications = NULL;
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
	struct rebx_node* current = rebx->forces;
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
    struct rebx_node* current = rebx->pre_timestep_operators;
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
    struct rebx_node* current = rebx->post_timestep_operators;
    const double dt = sim->dt_last_done;
    
    while(current != NULL){
        struct rebx_step* step = current->object;
        struct rebx_operator* operator = step->operator;
        operator->step(sim, operator, dt*step->dt_fraction);
        current = current->next;
    }
}

/**********************************************
 Adders for linked lists in extras structure
 *********************************************/

static struct rebx_node* rebx_create_node(struct rebx_extras* rebx){
    struct rebx_node* node = rebx_malloc(rebx, sizeof(*node));
    if (node == NULL){
        return NULL;
    }
    node->object = NULL;
    node->name = NULL;
    node->next = NULL;
    return node;
}

struct rebx_force* rebx_create_force(struct rebx_extras* const rebx, const char* name){
    struct rebx_force* force = rebx_malloc(rebx, sizeof(*force));
    if (force == NULL){
        return NULL;
    }
    force->ap = NULL;
    force->_sim = rebx->sim;
    force->name = rebx_malloc(rebx, strlen(name) + 1); // +1 for \0 at end
    if (force->name == NULL){
        return NULL;
    }
    else{
        strcpy(force->name, name);
    }
    
    // loop over all names. precompiler list
    if(strcmp(name, "gr") == 0){
        force->update_accelerations = rebx_gr;
        force->force_type = REBX_FORCE_VEL;
    }
    /*if (strcmp(name, "modify_orbits_forces") == 0){
        force->update_accelerations = rebx_modify_orbits_forces;
        force->force_type = REBX_FORCE_VEL;
    }
    
    else if (strcmp(name, "gr_full") == 0){
        force->update_accelerations = rebx_gr_full;
        force->force_type = REBX_FORCE_VEL;
    }
    else if (strcmp(name, "gravitational_harmonics") == 0){
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
        char str[100];
        sprintf(str, "REBOUNDx error: Force '%s' not found in rebx_create_force.\n", name);
        reb_error(rebx->sim, str);
        return NULL;
    }
    return force;
}

struct rebx_operator* rebx_create_operator(struct rebx_extras* const rebx, const char* name){
    struct rebx_operator* operator = rebx_malloc(rebx, sizeof(*operator));
    if (operator == NULL){
        return NULL;
    }
    operator->ap = NULL;
    operator->_sim = rebx->sim;
    operator->name = rebx_malloc(rebx, strlen(name) + 1); // +1 for \0 at end
    if (operator->name == NULL){
        return NULL;
    }
    else{
        strcpy(operator->name, name);
    }
    
    if (strcmp(name, "kepler") == 0){
        operator->step = rebx_kepler_step;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }
    else if (strcmp(name, "jump") == 0){
        operator->step = rebx_jump_step;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }
    else if (strcmp(name, "interaction") == 0){
        operator->step = rebx_interaction_step;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }
    else if (strcmp(name, "ias15") == 0){
        operator->step = rebx_ias15_step;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }
    /*
    else if (strcmp(name, "modify_orbits_direct") == 0){
        operator->step = rebx_modify_orbits_direct;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }
    else if (strcmp(name, "modify_mass") == 0){
        operator->step = rebx_modify_mass;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }*/
    
    /*else if (hash == reb_hash("track_min_distance")){
     operator->step = rebx_track_min_distance;
     }*/
    else{
        char str[100];
        sprintf(str, "REBOUNDx error: Operator '%s' not found in rebx_create_operator.\n", name);
        reb_error(rebx->sim, str);
        return NULL;
    }
    return operator;
}

int rebx_add_force(struct rebx_extras* rebx, struct rebx_force* force){
    if (force == NULL){
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
    node->name = force->name;
    rebx_add_node(&rebx->forces, node);
    
    return 1;
}

int rebx_add_operator_step(struct rebx_extras* rebx, struct rebx_operator* operator, const double dt_fraction, enum rebx_timing timing, char* name){
    struct rebx_step* step = rebx_malloc(rebx, sizeof(*step));
    if(step == NULL){
        return 0;
    }
    step->name = name;
    step->operator = operator;
    step->dt_fraction = dt_fraction;
    
    struct rebx_node* node = rebx_create_node(rebx);
    if (node == NULL){
        return 0;
    }
    node->object = step;
    node->name = step->name;
    if (timing == REBX_TIMING_PRE){
        rebx_add_node(&rebx->pre_timestep_operators, node);
        return 1;
    }
    if (timing == REBX_TIMING_POST){
        rebx_add_node(&rebx->post_timestep_operators, node);
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
        int success = rebx_add_operator_step(rebx, operator, dt_fraction, REBX_TIMING_POST, operator->name);
        return success;
    }
    
    switch(sim->integrator){
        case REB_INTEGRATOR_IAS15:
        // don't add pre-timestep b/c don't know what IAS will choose as dt
        {
            dt_fraction = 1.;
            int success = rebx_add_operator_step(rebx, operator, dt_fraction, REBX_TIMING_POST, operator->name);
            return success;
        }
        case REB_INTEGRATOR_WHFAST: // half step pre and post
        {
            dt_fraction = 1./2.;
            int success1 = rebx_add_operator_step(rebx, operator, dt_fraction, REBX_TIMING_PRE, operator->name);
            int success2 = rebx_add_operator_step(rebx, operator, dt_fraction, REBX_TIMING_POST, operator->name);
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

// put incompatible integrator warnings in add_force and add_operator
/*struct rebx_force* rebx_init_force(struct rebx_extras* rebx, const char* name){
    struct rebx_force* force = malloc(sizeof(*force));
    uint32_t hash = reb_hash(name);
    if (hash == reb_hash("modify_orbits_forces")){
        force->update_accelerations = rebx_modify_orbits_forces;
        force->force_type = REBX_EFFECT_FORCE_VEL;
    }
    else if(hash == reb_hash("gr")){
        force->update_accelerations = rebx_gr;
        force->force_type = REBX_EFFECT_FORCE_VEL;
    }*/
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
/*    else{
        char str[100];
        sprintf(str, "REBOUNDx error: Effect '%s' not found in rebx_init_force.\n", name);
        reb_error(rebx->sim, str);
        free(force);
        return;
    }
    
    force->effect = malloc(sizeof(*force->effect));
    force->effect->sim = rebx->sim;
    force->effect->ap = NULL;
    
    return force;
}*/
    
/*void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param){
    struct rebx_param_to_be_freed* newparam = malloc(sizeof(*newparam));
    newparam->param = param;

    newparam->next = rebx->params_to_be_freed;
    rebx->params_to_be_freed = newparam;
}*/

/*********************************************************************************
 Internal functions for dealing with parameters
 ********************************************************************************/

/*********************************************************************************
 User interface for parameters
 ********************************************************************************/

/*int rebx_remove_param(struct rebx_node** apptr, const char* const param_name){
    return rebx_remove_node(apptr, reb_hash(param_name));
}*/

/* Checks whether name is registered, and if so creates new param struct and adds it to apptr linked list.
   Returns pointer to new param, caller must allocate and set param->value.
 */
                                                   
static struct rebx_param* rebx_new_param(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* name){
    enum rebx_param_type type = rebx_get_type(rebx, name);
    if (type == REBX_TYPE_NONE){
        char str[300];
        sprintf(str, "REBOUNDx Error: Need to register parameter name '%s' before using it. See examples.\n", name);
        reb_error(rebx->sim, str);
        return NULL;
    }
    
    // Allocate and initialize new param struct
    struct rebx_param* param = rebx_malloc(rebx, sizeof(*param));
    if (param == NULL){
        return NULL;
    }
    param->type = type;
    param->name = rebx_malloc(rebx, strlen(name) + 1); // +1 for \0 at end
    if (param->name == NULL){
        return NULL;
    }
    else{
        strcpy(param->name, name);
    }
    
    // Wrap in node and add it to linked list
    struct rebx_node* node = rebx_create_node(rebx);
    if (node == NULL){
        return NULL;
    }
    node->object = param;
    node->name = param->name;
    rebx_add_node(apptr, node);
    
    return param;
}

int rebx_set_param_pointer(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, void* val){
    if (apptr == NULL){
        reb_error(rebx->sim, "REBOUNDx Error: Passed NULL apptr to rebx_add_param. See examples.\n");
        return 0;
    }
    // Check whether it already exists in linked list
    struct rebx_param* param = rebx_get_param_struct(rebx, *apptr, param_name);
    
    if(param == NULL){
        param = rebx_new_param(rebx, apptr, param_name);
        if (param == NULL){  // adding new param failed
            return 0;
        }
    }
    
    param->value = val;
    return 1;
}

int rebx_set_param_double(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, double val){
    if (apptr == NULL){
        reb_error(rebx->sim, "REBOUNDx Error: Passed NULL apptr to rebx_add_param. See examples.\n");
        return 0;
    }
    // Check whether it already exists in linked list
    struct rebx_param* param = rebx_get_param_struct(rebx, *apptr, param_name);
    
    if(param == NULL){      // Make new param and allocate memory
        param = rebx_new_param(rebx, apptr, param_name);
        if (param == NULL){ // adding new param failed
            return 0;
        }
        param->value = rebx_malloc(rebx, sizeof(double));
        double* valptr = param->value;
        *valptr = val;
    }
    else{                   // Update existing. Don't allocate
        double* valptr = param->value;
        *valptr = val;
    }
    
    return 1;
}

int rebx_set_param_int(struct rebx_extras* const rebx, struct rebx_node** apptr, const char* const param_name, int val){
    if (apptr == NULL){
        reb_error(rebx->sim, "REBOUNDx Error: Passed NULL apptr to rebx_add_param. See examples.\n");
        return 0;
    }
    // Check whether it already exists in linked list
    struct rebx_param* param = rebx_get_param_struct(rebx, *apptr, param_name);
    
    if(param == NULL){      // Make new param and allocate memory
        param = rebx_new_param(rebx, apptr, param_name);
        if (param == NULL){ // adding new param failed
            return 0;
        }
        param->value = rebx_malloc(rebx, sizeof(int));
        int* valptr = param->value;
        *valptr = val;
    }
    else{                   // Update existing. Don't allocate
        int* valptr = param->value;
        *valptr = val;
    }
    
    return 1;
}

struct rebx_param* rebx_get_param_struct(struct rebx_extras* rebx, struct rebx_node* ap, const char* const param_name){
    struct rebx_node* node = rebx_get_node(ap, param_name);
    if (node == NULL){ // don't want warnings for optional parameters so don't reb_error
        return NULL;
    }
    else{
        return node->object;
    }
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

/*void* rebx_get_param(struct rebx_extras* const rebx, struct rebx_node* ap, const char* const param_name){
    struct rebx_param* param = rebx_get_param_struct(ap, param_name);
    if (param == NULL){
        return NULL;
    }
    return param->value;
}

void* rebx_get_param_old(struct rebx_extras* const rebx, struct rebx_node* ap, const char* const param_name, enum rebx_param_type param_type){
    struct rebx_param* param = rebx_get_param_struct(ap, param_name);
    if (param == NULL){
        return NULL;
    }
    if (param->param_type != param_type){
        char str[300];
        sprintf(str, "REBOUNDx Error: Parameter '%s' passed to rebx_get_param_check was found but was of wrong type.  See documentation for your particular effect.  In python, you might need to add a dot at the end of the number when assigning a parameter that REBOUNDx expects as a float.\n", param_name);
        reb_error(rebx->sim, str);
        return NULL;
    }

    return param->value;
}*/

struct rebx_effect* rebx_get_effect(struct rebx_extras* const rebx, const char* const effect_name){
    struct rebx_effect* current;// = rebx->forces->effect; // stopgap fix!
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

enum rebx_param_type rebx_get_type(struct rebx_extras* rebx, const char* name){
    struct rebx_node* node = rebx_get_node(rebx->registered_params, name);
    
    if (node == NULL){ // param not found
        return REBX_TYPE_NONE;
    }
    
    enum rebx_param_type* type = node->object;
    
    return *type;
}

void rebx_register_param(struct rebx_extras* const rebx, const char* name, enum rebx_param_type type){
    
    struct rebx_node* node = rebx_get_node(rebx->registered_params, name);
    
    if (node != NULL){
        char str[300];
        sprintf(str, "REBOUNDx Error: Parameter name '%s' already in registered list. Cannot add duplicates.\n", name);
        reb_error(rebx->sim, str);
        return;
    }
    else{
        node = rebx_create_node(rebx);
    }
    
    if (node == NULL){
        return;
    }
    
    node->object = rebx_malloc(rebx, sizeof(type));
    memcpy(node->object, &type, sizeof(type));
    node->next = NULL;
    node->name = rebx_malloc(rebx, strlen(name) + 1); // +1 for \0 at end
    if (node->name == NULL){
        return;
    }
    else{
        strcpy(node->name, name);
    }
    
    rebx_add_node(&rebx->registered_params, node);
    return;
}

void* rebx_malloc(struct rebx_extras* const rebx, size_t memsize){
    void* ptr = malloc(memsize);
    if (ptr == NULL && memsize>0){
        reb_error(rebx->sim, "REBOUNDx Error: Could not allocate memory.\n");
        return NULL;
    }
    
    struct rebx_node* node = malloc(sizeof(*node));
    if (node == NULL){
        reb_error(rebx->sim, "REBOUNDx Error: Could not allocate memory.\n");
        return NULL;
    }
    node->object = ptr;
    node->next = NULL;
    
    rebx_add_node(&rebx->allocated_pointers, node);
    return ptr;
}

void rebx_free(struct rebx_extras* rebx){
    rebx_free_pointers(rebx);
    free(rebx);
}

void rebx_free_pointers(struct rebx_extras* rebx){
    struct rebx_node* current = rebx->allocated_pointers;
    struct rebx_node* next;
    while (current != NULL){
        next = current->next;
        free(current->object);
        free(current);
        current = next;
    }
}
