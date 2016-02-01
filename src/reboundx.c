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
const char* rebx_version_str = "2.4.1";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.

/* Main routines called each timestep. */

void rebx_forces(struct reb_simulation* sim){
    struct rebx_extras* rebx = sim->extras;
    for(int i=0; i<rebx->Nforces; i++){
        rebx->forces[i](sim);
    }
}

void rebx_ptm(struct reb_simulation* sim){
    struct rebx_extras* rebx = sim->extras;
    for(int i=0; i<rebx->Nptm; i++){
        rebx->ptm[i](sim);
    }
}

/* Initialization routine. */


void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx){
    sim->extras = rebx;
    rebx->sim = sim;

    rebx->params_to_be_freed = NULL;
    rebx->ptm = NULL;
    rebx->forces = NULL;

    rebx->modify_orbits_forces.p = 0;
    rebx->modify_orbits_forces.coordinates = JACOBI;
    rebx->modify_orbits_direct.p = 0;
    rebx->modify_orbits_direct.coordinates = JACOBI;
    rebx->gr.c = C_DEFAULT; 
    rebx->radiation_forces.c = C_DEFAULT;
    rebx->radiation_forces.source = NULL;

    rebx->Nptm = 0;
    rebx->Nforces = 0;
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

void rebx_free_pointers(struct rebx_extras* rebx){
    rebx_free_params(rebx);
    free(rebx->forces);
    free(rebx->ptm);
}


/* Internal utility functions. */

void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param){
    struct rebx_param_to_be_freed* newparam = malloc(sizeof(*newparam));
    newparam->param = param;

    newparam->next = rebx->params_to_be_freed;
    rebx->params_to_be_freed = newparam;
}

void* rebx_search_param(const struct reb_particle* p, enum REBX_PARAMS param){
    struct rebx_param* current = p->ap;
    while(current != NULL){
        if(current->param_type == param){
            return current->paramPtr;
        }
        current = current->next;
    }
    return NULL;
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

void rebx_add_param_orb_tau(struct reb_particle* p){
    struct rebx_param* newparam = malloc(sizeof(*newparam));
    newparam->paramPtr = malloc(sizeof(struct rebx_orb_tau));
    struct rebx_orb_tau orb_tau = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY}; // set all timescales to infinity (i.e. no effect)
    *(struct rebx_orb_tau*) newparam->paramPtr = orb_tau;
    newparam->param_type = (enum REBX_PARAMS)ORB_TAU;

    newparam->next = p->ap;
    p->ap = newparam;

    rebx_add_param_to_be_freed(p->sim->extras, newparam);
}

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

void rebx_add_gr(struct rebx_extras* rebx, double c){
    struct reb_simulation* sim = rebx->sim;
    sim->additional_forces = rebx_forces;
    sim->force_is_velocity_dependent = 1;

    rebx->Nforces++;
    rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
    rebx->forces[rebx->Nforces-1] = rebx_gr;
    rebx->gr.c = c;
}

void rebx_add_gr_full(struct rebx_extras* rebx, double c){
    struct reb_simulation* sim = rebx->sim;
    sim->additional_forces = rebx_forces;
    sim->force_is_velocity_dependent = 1;

    rebx->Nforces++;
    rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
    rebx->forces[rebx->Nforces-1] = rebx_gr_full;
    rebx->gr.c = c;
}

void rebx_add_gr_potential(struct rebx_extras* rebx, double c){
    struct reb_simulation* sim = rebx->sim;
    sim->additional_forces = rebx_forces;

    rebx->Nforces++;
    rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
    rebx->forces[rebx->Nforces-1] = rebx_gr_potential;
    rebx->gr.c = c;
}

void rebx_add_modify_orbits_forces(struct rebx_extras* rebx){
    struct reb_simulation* sim = rebx->sim;
    sim->additional_forces = rebx_forces;
    sim->force_is_velocity_dependent = 1;

    rebx->Nforces++;
    rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
    rebx->forces[rebx->Nforces-1] = rebx_modify_orbits_forces;
}

void rebx_add_modify_orbits_direct(struct rebx_extras* rebx){
    struct reb_simulation* sim = rebx->sim;
    sim->post_timestep_modifications = rebx_ptm;

    rebx->Nptm++;
    rebx->ptm = realloc(rebx->ptm, sizeof(*rebx->ptm)*rebx->Nptm);
    rebx->ptm[rebx->Nptm-1] = rebx_modify_orbits_direct;
}

void rebx_add_radiation_forces(struct rebx_extras* rebx, struct reb_particle* source, double c){
    struct reb_simulation* sim = rebx->sim;
    sim->additional_forces = rebx_forces;
    sim->force_is_velocity_dependent = 1;

    rebx->Nforces++;
    rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
    rebx->forces[rebx->Nforces-1] = rebx_radiation_forces;
    rebx->radiation_forces.source = source;
    rebx->radiation_forces.c = c;
}

void rebx_add_custom_ptm(struct rebx_extras* rebx, void (*custom_ptm)(struct reb_simulation* const r) ){
    struct reb_simulation* sim = rebx->sim;
    sim->post_timestep_modifications = rebx_ptm;

    rebx->Nptm++;
    rebx->ptm = realloc(rebx->ptm, sizeof(*rebx->ptm)*rebx->Nptm);
    rebx->ptm[rebx->Nptm-1] = custom_ptm;
}

void rebx_add_custom_forces(struct rebx_extras* rebx, void (*custom_forces)(struct reb_simulation* const r), int force_is_velocity_dependent){
    struct reb_simulation* sim = rebx->sim;
    sim->additional_forces = rebx_forces;
    if (force_is_velocity_dependent){
        sim->force_is_velocity_dependent = 1;
    }

    rebx->Nforces++;
    rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
    rebx->forces[rebx->Nforces-1] = custom_forces;
}

/****************************************
Functions for getting and setting particle parameters
 *****************************************/

// Getter setter landmark for add_param.py
void rebx_set_tau_a(struct reb_particle* p, double value){
    struct rebx_orb_tau* orb_tau = rebx_search_param(p, ORB_TAU);
    if(orb_tau == NULL){
        rebx_add_param_orb_tau(p);
        orb_tau = rebx_search_param(p, ORB_TAU);
    }
    orb_tau->tau_a = value;
}

void rebx_set_tau_e(struct reb_particle* p, double value){
    struct rebx_orb_tau* orb_tau = rebx_search_param(p, ORB_TAU);
    if(orb_tau == NULL){
        rebx_add_param_orb_tau(p);
        orb_tau = rebx_search_param(p, ORB_TAU);
    }
    orb_tau->tau_e = value;
}

void rebx_set_tau_inc(struct reb_particle* p, double value){
    struct rebx_orb_tau* orb_tau = rebx_search_param(p, ORB_TAU);
    if(orb_tau == NULL){
        rebx_add_param_orb_tau(p);
        orb_tau = rebx_search_param(p, ORB_TAU);
    }
    orb_tau->tau_inc = value;
}

void rebx_set_tau_omega(struct reb_particle* p, double value){
    struct rebx_orb_tau* orb_tau = rebx_search_param(p, ORB_TAU);
    if(orb_tau == NULL){
        rebx_add_param_orb_tau(p);
        orb_tau = rebx_search_param(p, ORB_TAU);
    }
    orb_tau->tau_omega = value;
}

void rebx_set_tau_Omega(struct reb_particle* p, double value){
    struct rebx_orb_tau* orb_tau = rebx_search_param(p, ORB_TAU);
    if(orb_tau == NULL){
        rebx_add_param_orb_tau(p);
        orb_tau = rebx_search_param(p, ORB_TAU);
    }
    orb_tau->tau_Omega = value;
}

double rebx_get_tau_a(struct reb_particle* p){
    struct rebx_orb_tau* orb_tau = rebx_search_param(p, ORB_TAU);
    if(orb_tau == NULL){
        return INFINITY;
    }
    else{
        return orb_tau->tau_a;
    }
}

double rebx_get_tau_e(struct reb_particle* p){
    struct rebx_orb_tau* orb_tau = rebx_search_param(p, ORB_TAU);
    if(orb_tau == NULL){
        return INFINITY;
    }
    else{
        return orb_tau->tau_e;
    }
}

double rebx_get_tau_inc(struct reb_particle* p){
    struct rebx_orb_tau* orb_tau = rebx_search_param(p, ORB_TAU);
    if(orb_tau == NULL){
        return INFINITY;
    }
    else{
        return orb_tau->tau_inc;
    }
}

double rebx_get_tau_omega(struct reb_particle* p){
    struct rebx_orb_tau* orb_tau = rebx_search_param(p, ORB_TAU);
    if(orb_tau == NULL){
        return INFINITY;    
    }
    else{
        return orb_tau->tau_omega;
    }
}

double rebx_get_tau_Omega(struct reb_particle* p){
    struct rebx_orb_tau* orb_tau = rebx_search_param(p, ORB_TAU);
    if(orb_tau == NULL){
        return INFINITY;
    }
    else{
        return orb_tau->tau_Omega;
    }
}

void rebx_set_beta(struct reb_particle* p, double value){
    double* betaPtr = rebx_search_param(p, RAD_BETA);
    if(betaPtr == NULL){
        rebx_add_param_double(p, RAD_BETA, value);
    }
    else{
        *betaPtr = value;
    }
}

double rebx_get_beta(struct reb_particle* p){
    double* betaPtr = rebx_search_param(p, RAD_BETA);
    if(betaPtr == NULL){
        return 0.;
    }
    else{
        return *betaPtr;
    }
}

/****************************************
Convenience Functions (include modification in function name in some form)
 *****************************************/

double rebx_rad_calc_beta(struct rebx_extras* rebx, double particle_radius, double density, double Q_pr, double L){
    double mu = rebx->sim->G*rebx->radiation_forces.source->m;
    const double c = rebx->radiation_forces.c;
    return 3.*L*Q_pr/(16.*M_PI*mu*c*density*particle_radius);   
}
double rebx_rad_calc_particle_radius(struct rebx_extras* rebx, double beta, double density, double Q_pr, double L){
    const double mu = rebx->sim->G*rebx->radiation_forces.source->m;
    const double c = rebx->radiation_forces.c;
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
