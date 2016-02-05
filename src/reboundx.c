/**
 * @file    reboundx.c
 * @brief   Initialization and getter/setter routines.
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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include "reboundx.h"

const char* rebx_build_str = __DATE__ " " __TIME__; // Date and time build string. 
const char* rebx_version_str = "2.5.0";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.

/* Main routines called each timestep. */

void rebx_forces(struct reb_simulation* sim){
	struct rebx_extras* rebx = sim->extras;
	struct rebx_effect* current = rebx->forces;
	while(current != NULL){
		current->functionPtr(sim, current);
		current = current->next;
	}
}

void rebx_post_timestep_modifications(struct reb_simulation* sim){
	struct rebx_extras* rebx = sim->extras;
	struct rebx_effect* current = rebx->post_timestep_modifications;
	while(current != NULL){
		current->functionPtr(sim, current);
		current = current->next;
	}
}

/* Initialization routine. */

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx){
    sim->extras = rebx;
    rebx->sim = sim;

	rebx->params_to_be_freed = NULL;
	rebx->post_timestep_modifications = NULL;
	rebx->forces = NULL;
}

/* Garbage collection routines. */

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

void rebx_free_effects(struct rebx_effect* effects){
    struct rebx_effect* current = effects;
    struct rebx_effect* temp_next;

    while(current != NULL){
        temp_next = current->next;
        free(current->paramsPtr);
        free(current);
        current = temp_next;
    }
}

void rebx_free_pointers(struct rebx_extras* rebx){
    rebx_free_params(rebx);
    rebx_free_effects(rebx->forces);
    rebx_free_effects(rebx->post_timestep_modifications);
}

/* Internal utility functions. */

void rebx_add_force(struct rebx_extras* rebx, enum REBX_EFFECTS effect_type, void* paramsPtr, void (*functionPtr) (struct reb_simulation* sim, struct rebx_effect* effect)){
    struct rebx_effect* effect = malloc(sizeof(*effect));
    effect->paramsPtr = paramsPtr;
    effect->functionPtr = functionPtr;
    effect->effect_type = effect_type;
    effect->next = rebx->forces;
    rebx->forces = effect;
}

void rebx_add_post_timestep_modification(struct rebx_extras* rebx, enum REBX_EFFECTS effect_type, void* paramsPtr, void (*functionPtr) (struct reb_simulation* sim, struct rebx_effect* effect)){
    struct rebx_effect* effect = malloc(sizeof(*effect));
    effect->paramsPtr = paramsPtr;
    effect->functionPtr = functionPtr;
    effect->effect_type = effect_type;
    effect->next = rebx->post_timestep_modifications;
    rebx->post_timestep_modifications = effect;
}

void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param){
    struct rebx_param_to_be_freed* newparam = malloc(sizeof(*newparam));
    newparam->param = param;

    newparam->next = rebx->params_to_be_freed;
    rebx->params_to_be_freed = newparam;
}

void* rebx_get_param(const struct reb_particle* p, enum REBX_PARAMS param){
    struct rebx_param* current = p->ap;
    while(current != NULL){
        if(current->param_type == param){
            return current->paramPtr;
        }
        current = current->next;
    }
    return NULL;
}

struct rebx_effect* rebx_get_effect_in(struct rebx_effect* effects, enum REBX_EFFECTS effect_type){
    struct rebx_effect* current = effects;
    while(current != NULL){
        if(current->effect_type == effect_type){
            return current;
        }
        current = current->next;
    }
    return NULL;
}

struct rebx_effect* rebx_get_effect(struct rebx_extras* rebx, enum REBX_EFFECTS effect_type){
    /* Effects can't be in both forces and post_timestep_modifications, so simply check both*/
    struct rebx_effect* effect = rebx_get_effect_in(rebx->forces, effect_type);
    if(effect == NULL){
        effect = rebx_get_effect_in(rebx->post_timestep_modifications, effect_type);
    }
    return effect;
}


void* rebx_get_effect_params_in(struct rebx_effect* effects, enum REBX_EFFECTS effect_type){
    struct rebx_effect* effect = rebx_get_effect_in(effects, effect_type);
    if(effect == NULL){
        return NULL;
    }
    else{
        return effect->paramsPtr;
    }
}

/* Internal parameter adders (need a different one for each REBX_PARAM type). */
/* Generic adder for params that are a single double value */
void rebx_add_param_double(struct reb_particle* p, enum REBX_PARAMS param_type, double value){
    struct rebx_param* newparam = malloc(sizeof(*newparam));
    newparam->paramPtr = malloc(sizeof(double));
    *(double*) newparam->paramPtr = value;
    newparam->param_type = param_type;

    newparam->next = p->ap;
    p->ap = newparam;

    rebx_add_param_to_be_freed(p->sim->extras, newparam);
}   

/*void rebx_add_param_orb_tau(struct reb_particle* p){
    struct rebx_param* newparam = malloc(sizeof(*newparam));
    newparam->paramPtr = malloc(sizeof(struct rebx_orb_tau));
    struct rebx_orb_tau orb_tau = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY}; // set all timescales to infinity (i.e. no effect)
    *(struct rebx_orb_tau*) newparam->paramPtr = orb_tau;
    newparam->param_type = (enum REBX_PARAMS)ORB_TAU;

    newparam->next = p->ap;
    p->ap = newparam;

    rebx_add_param_to_be_freed(p->sim->extras, newparam);
}*/

/****************************************
User API
*****************************************/

struct rebx_extras* rebx_init(struct reb_simulation* sim){
    struct rebx_extras* rebx = malloc(sizeof(*rebx));
    rebx_initialize(sim, rebx);
    return rebx;
}

void rebx_free(struct rebx_extras* rebx){
    rebx_free_pointers(rebx);
    free(rebx);
}

/****************************************
Functions for specific REBOUNDx effects
 *****************************************/

void* rebx_get_effect_params(struct rebx_extras* rebx, enum REBX_EFFECTS effect_type){
    /* Effects can't be in both forces and post_timestep_modifications, so simply check both*/
    void* params = rebx_get_effect_params_in(rebx->forces, effect_type);
    if(params == NULL){
        params = rebx_get_effect_params_in(rebx->post_timestep_modifications, effect_type);
    }
    return params;
}

void rebx_add_gr(struct rebx_extras* rebx, struct reb_particle* source, double c){
	struct reb_simulation* sim = rebx->sim;
	sim->additional_forces = rebx_forces;
	
	struct rebx_params_gr* gr_params = malloc(sizeof(*gr_params));
	gr_params->c = c;
    if(source == NULL){
        gr_params->source_index = 0;
    }
    else{
        gr_params->source_index = reb_get_particle_index(source);
    }
	
    sim->force_is_velocity_dependent = 1;
    rebx_add_force(rebx, GR, gr_params, rebx_gr);
}

void rebx_add_gr_full(struct rebx_extras* rebx, double c){
	struct reb_simulation* sim = rebx->sim;
	sim->additional_forces = rebx_forces;
	
	struct rebx_params_gr* gr_params = malloc(sizeof(*gr_params));
	gr_params->c = c;
	
    sim->force_is_velocity_dependent = 1;
    rebx_add_force(rebx, GR_FULL, gr_params, rebx_gr_full);
}

void rebx_add_gr_potential(struct rebx_extras* rebx, struct reb_particle* source, double c){
	struct reb_simulation* sim = rebx->sim;
	sim->additional_forces = rebx_forces;
	
	struct rebx_params_gr* gr_params = malloc(sizeof(*gr_params));
	gr_params->c = c;
    if(source == NULL){
        gr_params->source_index = 0;
    }
    else{
        gr_params->source_index = reb_get_particle_index(source);
    }
	
    rebx_add_force(rebx, GR_POTENTIAL, gr_params, rebx_gr_potential);
}

void rebx_add_modify_orbits_forces(struct rebx_extras* rebx){
	struct reb_simulation* sim = rebx->sim;
	sim->additional_forces = rebx_forces;
	
	struct rebx_params_modify_orbits* mod_params = malloc(sizeof(*mod_params));
	mod_params->p = 0;
    mod_params->coordinates = JACOBI;
	
    sim->force_is_velocity_dependent = 1;
    rebx_add_force(rebx, MODIFY_ORBITS_FORCES, mod_params, rebx_modify_orbits_forces);
}

void rebx_add_modify_orbits_direct(struct rebx_extras* rebx){
	struct reb_simulation* sim = rebx->sim;
	sim->post_timestep_modifications = rebx_post_timestep_modifications;
	
	struct rebx_params_modify_orbits* mod_params = malloc(sizeof(*mod_params));
	mod_params->p = 0;
    mod_params->coordinates = JACOBI;
	
    rebx_add_post_timestep_modification(rebx, MODIFY_ORBITS_DIRECT, mod_params, rebx_modify_orbits_direct);
}

void rebx_add_radiation_forces(struct rebx_extras* rebx, struct reb_particle* source, double c){
	struct reb_simulation* sim = rebx->sim;
	sim->additional_forces = rebx_forces;
	
	struct rebx_params_radiation_forces* rad_params = malloc(sizeof(*rad_params));
	rad_params->c = c;
    if(source == NULL){
        rad_params->source_index = 0;
    }
    else{
        rad_params->source_index = reb_get_particle_index(source);
    }
	
    sim->force_is_velocity_dependent = 1;
    rebx_add_force(rebx, RADIATION_FORCES, rad_params, rebx_radiation_forces);
}

void rebx_add_custom_post_timestep_modification(struct rebx_extras* rebx, void (*custom_post_timestep_modification)(struct reb_simulation* const sim, struct rebx_effect* custom_effect), void* custom_params){
    struct reb_simulation* sim = rebx->sim;
    sim->post_timestep_modifications = rebx_post_timestep_modifications;

    rebx_add_post_timestep_modification(rebx, CUSTOM_POST_TIMESTEP_MODIFICATION, custom_params, custom_post_timestep_modification);
}

void rebx_add_custom_force(struct rebx_extras* rebx, void (*custom_force)(struct reb_simulation* const sim, struct rebx_effect* custom_effect), int force_is_velocity_dependent, void* custom_params){
    struct reb_simulation* sim = rebx->sim;
    sim->additional_forces = rebx_forces;
    if (force_is_velocity_dependent){
        sim->force_is_velocity_dependent = 1;
    }
    rebx_add_force(rebx, CUSTOM_POST_TIMESTEP_MODIFICATION, custom_params, custom_force);
}

/****************************************************
Functions for getting and setting particle parameters
 ****************************************************/

// Getter setter landmark for add_param.py
void rebx_set_beta(struct reb_particle* p, double value){
    double* betaPtr = rebx_get_param(p, RAD_BETA);
    if(betaPtr == NULL){
        rebx_add_param_double(p, RAD_BETA, value);
    }
    else{
        *betaPtr = value;
    }
}

double rebx_get_beta(struct reb_particle* p){
    double* betaPtr = rebx_get_param(p, RAD_BETA);
    if(betaPtr == NULL){
        return 0.;
    }
    else{
        return *betaPtr;
    }
}

void rebx_set_tau_omega(struct reb_particle* p, double value){
    double* tau_omegaPtr = rebx_get_param(p, TAU_LITTLE_OMEGA);
    if(tau_omegaPtr == NULL){
        rebx_add_param_double(p, TAU_LITTLE_OMEGA, value);
    }
    else{
        *tau_omegaPtr = value;
    }
}

double rebx_get_tau_omega(struct reb_particle* p){
    double* tau_omegaPtr = rebx_get_param(p, TAU_LITTLE_OMEGA);
    if(tau_omegaPtr == NULL){
        return INFINITY;
    }
    else{
        return *tau_omegaPtr;
    }
}

void rebx_set_tau_Omega(struct reb_particle* p, double value){
    double* tau_OmegaPtr = rebx_get_param(p, TAU_BIG_OMEGA);
    if(tau_OmegaPtr == NULL){
        rebx_add_param_double(p, TAU_BIG_OMEGA, value);
    }
    else{
        *tau_OmegaPtr = value;
    }
}

double rebx_get_tau_Omega(struct reb_particle* p){
    double* tau_OmegaPtr = rebx_get_param(p, TAU_BIG_OMEGA);
    if(tau_OmegaPtr == NULL){
        return INFINITY;
    }
    else{
        return *tau_OmegaPtr;
    }
}

void rebx_set_tau_inc(struct reb_particle* p, double value){
    double* tau_incPtr = rebx_get_param(p, TAU_INC);
    if(tau_incPtr == NULL){
        rebx_add_param_double(p, TAU_INC, value);
    }
    else{
        *tau_incPtr = value;
    }
}

double rebx_get_tau_inc(struct reb_particle* p){
    double* tau_incPtr = rebx_get_param(p, TAU_INC);
    if(tau_incPtr == NULL){
        return INFINITY;
    }
    else{
        return *tau_incPtr;
    }
}

void rebx_set_tau_e(struct reb_particle* p, double value){
    double* tau_ePtr = rebx_get_param(p, TAU_E);
    if(tau_ePtr == NULL){
        rebx_add_param_double(p, TAU_E, value);
    }
    else{
        *tau_ePtr = value;
    }
}

double rebx_get_tau_e(struct reb_particle* p){
    double* tau_ePtr = rebx_get_param(p, TAU_E);
    if(tau_ePtr == NULL){
        return INFINITY;
    }
    else{
        return *tau_ePtr;
    }
}

void rebx_set_tau_a(struct reb_particle* p, double value){
    double* tau_aPtr = rebx_get_param(p, TAU_A);
    if(tau_aPtr == NULL){
        rebx_add_param_double(p, TAU_A, value);
    }
    else{
        *tau_aPtr = value;
    }
}

double rebx_get_tau_a(struct reb_particle* p){
    double* tau_aPtr = rebx_get_param(p, TAU_A);
    if(tau_aPtr == NULL){
        return INFINITY;
    }
    else{
        return *tau_aPtr;
    }
}

/****************************************
Convenience Functions (include modification in function name in some form)
 *****************************************/

double rebx_rad_calc_beta(struct rebx_extras* rebx, double particle_radius, double density, double Q_pr, double L){
    struct rebx_params_radiation_forces* params = rebx_get_effect_params_in(rebx->forces, RADIATION_FORCES);
    double mu = rebx->sim->G*rebx->sim->particles[params->source_index].m;
    const double c = params->c;
    return 3.*L*Q_pr/(16.*M_PI*mu*c*density*particle_radius);   
}
double rebx_rad_calc_particle_radius(struct rebx_extras* rebx, double beta, double density, double Q_pr, double L){
    struct rebx_params_radiation_forces* params = rebx_get_effect_params_in(rebx->forces, RADIATION_FORCES);
    double mu = rebx->sim->G*rebx->sim->particles[params->source_index].m;
    const double c = params->c;
    return 3.*L*Q_pr/(16.*M_PI*mu*c*density*beta);
}

/* Function to test whether REBOUNDx can load librebound.so and call REBOUND functions. */

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
