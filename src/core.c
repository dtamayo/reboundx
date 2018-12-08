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
    //rebx_free_params(rebx);
    rebx_free_effects(rebx);
    free(rebx);
}

/*void rebx_free_params(struct rebx_extras* rebx){
    struct rebx_param__to_be_freed* current = rebx->params_to_be_freed;
    struct rebx_param_to_be_freed* temp_next;
    while(current != NULL){
        temp_next = current->next;
        free(current->param->contents);
        free(current->param);
        free(current);
        current = temp_next;
    }
}*/

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
		force->update_accelerations(sim, force->effect, sim->particles, N);
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
        struct rebx_operator* operator = current->object;
        operator->step(sim, operator->effect, dt*operator->dtfactor);
        current = current->next;
    }
}

void rebx_post_timestep_modifications(struct reb_simulation* sim){
    struct rebx_extras* rebx = sim->extras;
    struct rebx_node* current = rebx->post_timestep_operators;
    const double dt = sim->dt_last_done;
    
    while(current != NULL){
        struct rebx_operator* operator = current->object;
        operator->step(sim, operator->effect, dt*operator->dtfactor);
        current = current->next;
    }
}

/**********************************************
 Adders for linked lists in extras structure
 *********************************************/

struct rebx_effect* rebx_create_effect(struct reb_simulation* const sim){
    struct rebx_effect* effect = rebx_malloc(sim, sizeof(*effect));
    if (effect == NULL){
        return NULL;
    }
    
    effect->ap = NULL;
    effect->sim = sim;
    return effect;
}

static struct rebx_force* rebx_create_force_silent(struct reb_simulation* const sim, const char* name, struct rebx_effect* effect){
    struct rebx_force* force = rebx_malloc(sim, sizeof(*force));
    if (force == NULL){
        return NULL;
    }
    force->effect = effect;
    if (strcmp(name, "modify_orbits_forces") == 0){
        force->update_accelerations = rebx_modify_orbits_forces;
        force->force_type = REBX_FORCE_VEL;
    }
    else if(strcmp(name, "gr") == 0){
        force->update_accelerations = rebx_gr;
        force->force_type = REBX_FORCE_VEL;
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
        free(force);
        return NULL;
    }
    return force;
}

struct rebx_force* rebx_create_force(struct reb_simulation* const sim, const char* name, struct rebx_effect* effect){
    struct rebx_force* force = rebx_create_force_silent(sim, name, effect);
    if (!force){
        char str[100];
        sprintf(str, "REBOUNDx error: Force '%s' not found in rebx_create_force.\n", name);
        reb_error(sim, str);
        return NULL;
    }
    return force;
}

static struct rebx_operator* rebx_create_operator_silent(struct reb_simulation* const sim, const char* name, const double dtfactor, struct rebx_effect* effect){
    struct rebx_operator* operator = rebx_malloc(sim, sizeof(*operator));
    if (operator == NULL){
        return NULL;
    }
    operator->effect = effect;
    operator->dtfactor = dtfactor;
    if (strcmp(name, "modify_orbits_direct") == 0){
        operator->step = rebx_modify_orbits_direct;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }
    else if (strcmp(name, "modify_mass") == 0){
        operator->step = rebx_modify_mass;
        operator->operator_type = REBX_OPERATOR_UPDATER;
    }
    /*else if (hash == reb_hash("track_min_distance")){
     operator->step = rebx_track_min_distance;
     }*/
    else{
        free(operator);
        return NULL;
    }
    return operator;
}

struct rebx_operator* rebx_create_operator(struct reb_simulation* const sim, const char* name, const double dtfactor, struct rebx_effect* effect){
    struct rebx_operator* operator = rebx_create_operator_silent(sim, name, dtfactor, effect);
    if (!operator){
        char str[100];
        sprintf(str, "REBOUNDx error: Operator '%s' not found in rebx_create_operator.\n", name);
        reb_error(sim, str);
        return NULL;
    }
    return operator;
}

void rebx_add_force(struct rebx_extras* rebx, const char* name, struct rebx_effect* effect){
    struct rebx_force* force = rebx_create_force(rebx->sim, name, effect);
    if (force == NULL){
        return;
    }
    if (force->force_type == REBX_FORCE_VEL){
        rebx->sim->force_is_velocity_dependent = 1;
    }
    rebx_add_node(rebx->sim, &rebx->forces, force, name);
}

void rebx_add_pre_timestep_operator(struct rebx_extras* rebx, const char* name, const double dtfactor, struct rebx_effect* effect){
    struct rebx_operator* operator = rebx_create_operator(rebx->sim, name, dtfactor, effect);
    if(operator == NULL){
        return;
    }
    rebx_add_node(rebx->sim, &rebx->pre_timestep_operators, operator, name);
}

void rebx_add_post_timestep_operator(struct rebx_extras* rebx, const char* name, const double dtfactor, struct rebx_effect* effect){
    struct rebx_operator* operator = rebx_create_operator(rebx->sim, name, dtfactor, effect);
    if(operator == NULL){
        return;
    }
    rebx_add_node(rebx->sim, &rebx->post_timestep_operators, operator, name);
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

struct rebx_effect* rebx_add(struct rebx_extras* rebx, const char* name){
    struct reb_simulation* const sim = rebx->sim;
    struct rebx_effect* effect = rebx_create_effect(sim);
    struct rebx_force* force = rebx_create_force_silent(sim, name, effect); // name might be an operator. Don't warn if not found.
    if (force){
        if (force->force_type == REBX_FORCE_VEL){
            rebx->sim->force_is_velocity_dependent = 1;
        }
        
        if(sim->integrator == REB_INTEGRATOR_IAS15){
            rebx_add_node(sim, &rebx->forces, force, name);
        }
        else if (sim->integrator == REB_INTEGRATOR_WHFAST){
            // Could add as operator if vel dependent here
            rebx_add_node(sim, &rebx->forces, force, name);
        }
        else if (sim->integrator == REB_INTEGRATOR_MERCURIUS){
            rebx_add_node(sim, &rebx->forces, force, name);
        }
        return force->effect;
    }
    
    double dtfactor = 1.;
    struct rebx_operator* operator = rebx_create_operator_silent(sim, name, dtfactor, effect);
    if (operator){
        if (sim->integrator == REB_INTEGRATOR_IAS15){
            rebx_add_node(sim, &rebx->post_timestep_operators, operator, name);
        }
        else if (sim->integrator == REB_INTEGRATOR_WHFAST){
            operator->dtfactor = 1./2.;
            rebx_add_node(sim, &rebx->pre_timestep_operators, operator, name);
            rebx_add_node(sim, &rebx->post_timestep_operators, operator, name);
        }
        else if (sim->integrator == REB_INTEGRATOR_MERCURIUS){
            if (operator->operator_type == REBX_OPERATOR_UPDATER){
                reb_error(sim, "REBOUNDx Error: Operators that change particle states are not supported with Mercurius.\n");
                free(operator);
                free(effect);
                return NULL;
            }
            else if (operator->operator_type == REBX_OPERATOR_RECORDER){
                // Doesn't alter state. Add once after timestep.
                rebx_add_node(sim, &rebx->post_timestep_operators, operator, name);
            }
        }
        return operator->effect;
    }
    char str[100];
    sprintf(str, "REBOUNDx error: Effect '%s' not found in rebx_add.\n", name);
    reb_error(sim, str);
    free(effect);
    return NULL;
}

struct rebx_effect* rebx_add_custom_force(struct rebx_extras* rebx, const char* name, void (*custom_force)(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N), const int force_is_velocity_dependent){
    struct rebx_effect* effect = rebx_create_effect(rebx->sim);
    /*effect->force = custom_force;
    struct reb_simulation* sim = rebx->sim;
    if(force_is_velocity_dependent){
        sim->force_is_velocity_dependent = 1;
    }*/
    return effect;
}

struct rebx_effect* rebx_add_custom_operator(struct rebx_extras* rebx, const char* name, void (*custom_operator)(struct reb_simulation* const sim, struct rebx_effect* const effect, const double dt, enum rebx_timing timing)){
    struct rebx_effect* effect = rebx_create_effect(rebx->sim);
    //effect->operator = custom_operator;
    return effect;
}
    
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

void* rebx_alloc_param_value(struct reb_simulation* const sim, enum rebx_param_type param_type){
    void* value = NULL;
    switch(param_type){
        case REBX_TYPE_DOUBLE:
        {
            value = rebx_malloc(sim, sizeof(double));
            break;
        }
        case REBX_TYPE_INT:
        {
            value = rebx_malloc(sim, sizeof(int));
            break;
        }
        case REBX_TYPE_UINT32:
        {
            value = rebx_malloc(sim, sizeof(uint32_t));
            break;
        }
        case REBX_TYPE_ORBIT:
        {
            value = rebx_malloc(sim, sizeof(struct reb_orbit));
            break;
        }
        case REBX_TYPE_LONGLONG:
        {
            value = rebx_malloc(sim, sizeof(long long));
            break;
        }
        // Nested structures not allocated so user does not mistakenly access unallocated fields. Treat together
        case REBX_TYPE_EFFECT:
        case REBX_TYPE_OPERATOR:
        case REBX_TYPE_FORCE:
        {
            reb_error(sim, "REBOUNDx Error: Nested structures cannot be allocated in rebx_alloc_param.\n");
            return NULL;
        }
        default:
        {
            char str[300];
            sprintf(str, "REBOUNDx Error: Parameter type '%d' passed to rebx_alloc_param not supported.\n", param_type);
            reb_error(sim, str);
            return NULL;
        }
    }
    
    return value;
}

// Must allocate memory for value first. Necessary for complicated nested parameters
// create allocates and initializes
struct rebx_param* rebx_create_param(struct reb_simulation* const sim, void* value, enum rebx_param_type param_type){
    struct rebx_param* param = rebx_malloc(sim, sizeof(*param));
    if (param == NULL){
        return NULL;
    }
    
    param->param_type = param_type;
    param->python_type = -1; // not used by C
    param->value = value;
    return param;
}

// add functions always require a linked list head
// need to allocate from inside out in nested structs to make sure you have refs to all pointers
void* rebx_add_param(struct reb_simulation* const sim, struct rebx_node** apptr, const char* const param_name, enum rebx_param_type param_type){
    if (apptr == NULL){
        reb_error(sim, "REBOUNDx Error: Passed NULL apptr to rebx_add_param. See examples.\n");
        return NULL;
    }
    
    // Check it doesn't already exist in linked list
    struct rebx_node* node = rebx_get_node(*apptr, reb_hash(param_name));
    if (node != NULL){
        char str[300];
        sprintf(str, "REBOUNDx Error: Parameter '%s' passed to rebx_add_param already exists.\n", param_name);
        reb_error(sim, str);
        return NULL;
    }
    
    // Doesn't exist, allocate and add new one
    void* value = rebx_alloc_param_value(sim, param_type);
    if (value == NULL){
        return NULL;
    }
    
    struct rebx_param* param = rebx_create_param(sim, value, param_type);
    if (param == NULL){
        free(value);
        return NULL;
    }
    
    node = rebx_add_node(sim, apptr, param, param_name);
    if (node){
        return param->value;
    }
    else{
        free(value);
        free(param);
        return NULL;
    }
}

struct rebx_param* rebx_get_param_struct(struct rebx_node* ap, const char* const param_name){
    struct rebx_node* node = rebx_get_node(ap, reb_hash(param_name));
    if (node == NULL){
        return NULL;
    }
    struct rebx_param* param = node->object;
    return param;
}

void* rebx_get_param(struct rebx_node* ap, const char* const param_name){
    struct rebx_param* param = rebx_get_param_struct(ap, param_name);
    if (param == NULL){
        return NULL;
    }
    return param->value;
}

void* rebx_get_param_check(struct reb_simulation* sim, struct rebx_node* ap, const char* const param_name, enum rebx_param_type param_type){
    struct rebx_param* param = rebx_get_param_struct(ap, param_name);
    if (param == NULL){
        return NULL;
    }
    if (param->param_type != param_type){
        char str[300];
        sprintf(str, "REBOUNDx Error: Parameter '%s' passed to rebx_get_param_check was found but was of wrong type.  See documentation for your particular effect.  In python, you might need to add a dot at the end of the number when assigning a parameter that REBOUNDx expects as a float.\n", param_name);
        reb_error(sim, str);
        return NULL;
    }

    return param->value;
}

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

void* rebx_malloc(struct reb_simulation* const sim, size_t memsize){
    void* ptr = malloc(memsize);
    if (ptr == NULL && memsize>0){
        reb_error(sim, "REBOUNDx Error: Could not allocate memory.\n");
        return NULL;
    }
    return ptr;
}
