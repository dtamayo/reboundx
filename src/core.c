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
#include "core.h"
#include "rebound.h"

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

void rebx_add_effect(struct rebx_extras* rebx, const char* name){
    struct rebx_effect* effect = malloc(sizeof(*effect));
    effect->next = rebx->effects;
    rebx->effects = effect;
    
    struct reb_simulation* sim = rebx->sim;
    switch(reb_tools_hash(name)){
        case reb_tools_hash("gr"):
            sim->force_is_velocity_dependent = 1;
            effect.force = rebx_gr;
            break;
        default:
            fprintf(stderr, "Effect not found.\n");
            exit(1);
    }
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

void* rebx_get_ap_hash(struct rebx_param current, uint32_t hash){
    while(current != NULL){
        if(current->hash == hash){
            return current->paramPtr;
        }
        current = current->next;
    }
    return NULL;
}

/*********************************************************************************
 Getters and Setters for particle parameters (need new set for each variable type)
 ********************************************************************************/

void rebx_set_ap_double(struct rebx_extras* rebx, struct rebx_param* ap, const char* ap_name, double value){
    uint32_t hash = reb_tools_hash(ap_name);
    rebx_set_ap_double_hash(rebx, ap, hash, value);
}

void rebx_set_ap_double_hash(struct rebx_extras* rebx, struct rebx_param* ap, uint32_t hash, double value){
    double* ptr = rebx_get_param_hash(ap, hash);
    if(ptr == NULL){
        ptr = rebx_add_ap_double(rebx, ap, hash);   // add a new parameter if it doesn't exist
    *ptr = value;                                   // update existing value
}

double* rebx_add_ap_double(struct rebx_extras* rebx, struct rebx_param* ap, uint32_t hash){
    struct rebx_param* newparam = malloc(sizeof(*newparam));
    newparam->paramPtr = malloc(sizeof(double));
    newparam->hash = hash;

    newparam->next = ap;
    ap = newparam;
    rebx_add_param_to_be_freed(rebx, newparam);
    return newparam;
}   

double rebx_get_ap_double(struct rebx_param* ap, const char* ap_name){
    uint32_t hash = reb_tools_hash(ap_name);
    return rebx_get_ap_double_hash(p, hash);
}

double rebx_get_ap_double_hash(struct rebx_param* ap, uint32_t hash){
    double* ptr = rebx_get_param_hash(p, hash);
    if(ptr == NULL){
        return NAN;
    }else{
        return *ptr;
    }
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
