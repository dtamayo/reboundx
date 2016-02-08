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
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include "reboundx.h"

const char* rebx_build_str = __DATE__ " " __TIME__; // Date and time build string. 
const char* rebx_version_str = "2.5.0";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.

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

/*void* rebx_get_effect_params(struct rebx_extras* rebx, enum REBX_EFFECTS effect_type){
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
*/
/****************************************************
Functions for getting and setting particle parameters
 ****************************************************/

// Getter setter landmark for add_param.py
/*void rebx_set_beta(struct reb_particle* p, double value){
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
*/
/****************************************
Convenience Functions (include modification in function name in some form)
 *****************************************/

/*double rebx_rad_calc_beta(struct rebx_extras* rebx, double particle_radius, double density, double Q_pr, double L){
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
*/
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
