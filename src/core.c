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
const char* rebx_version_str = "2.8.4";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.


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
        free(current->paramsPtr);
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
        if(current->is_force == 1){
		    current->functionPtr(sim, current);
        }
		current = current->next;
	}
}

void rebx_post_timestep_modifications(struct reb_simulation* sim){
	struct rebx_extras* rebx = sim->extras;
	struct rebx_effect* current = rebx->effects;
	while(current != NULL){
        if(current->is_force == 0){
		    current->functionPtr(sim, current);
        }
		current = current->next;
	}
}

/**********************************************
 Adders for linked lists in extras structure
 *********************************************/

void rebx_add_force(struct rebx_extras* rebx, void* paramsPtr, const char* name, void (*functionPtr) (struct reb_simulation* sim, struct rebx_effect* effect), int force_is_velocity_dependent){
    struct reb_simulation* sim = rebx->sim;
    sim->additional_forces = rebx_forces;
    if (force_is_velocity_dependent){
        sim->force_is_velocity_dependent = 1;
    }
    
    struct rebx_effect* effect = malloc(sizeof(*effect));
    effect->paramsPtr = paramsPtr;
    effect->functionPtr = functionPtr;
    effect->effect_type = rebx_hash(name);
    effect->is_force = 1;
    effect->next = rebx->effects;
    rebx->effects = effect;
}

void rebx_add_post_timestep_modification(struct rebx_extras* rebx, void* paramsPtr, const char* name, void (*functionPtr) (struct reb_simulation* sim, struct rebx_effect* effect)){
    struct reb_simulation* sim = rebx->sim;
    sim->post_timestep_modifications = rebx_post_timestep_modifications;

    struct rebx_effect* effect = malloc(sizeof(*effect));
    effect->paramsPtr = paramsPtr;
    effect->functionPtr = functionPtr;
    effect->effect_type = rebx_hash(name);
    effect->is_force = 0;
    effect->next = rebx->effects;
    rebx->effects = effect;
}

void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param){
    struct rebx_param_to_be_freed* newparam = malloc(sizeof(*newparam));
    newparam->param = param;

    newparam->next = rebx->params_to_be_freed;
    rebx->params_to_be_freed = newparam;
}

/****************************************
 Hash function
 *****************************************/
#define ROT32(x, y) ((x << y) | (x >> (32 - y))) // avoid effort
uint32_t rebx_murmur3_32(const char *key, uint32_t len, uint32_t seed) {
    // Source: Wikipedia
    static const uint32_t c1 = 0xcc9e2d51;
    static const uint32_t c2 = 0x1b873593;
    static const uint32_t r1 = 15;
    static const uint32_t r2 = 13;
    static const uint32_t m = 5;
    static const uint32_t n = 0xe6546b64;

    uint32_t hash = seed;

    const int nblocks = len / 4;
    const uint32_t *blocks = (const uint32_t *) key;
    int i;
    uint32_t k;
    for (i = 0; i < nblocks; i++) {
        k = blocks[i];
        k *= c1;
        k = ROT32(k, r1);
        k *= c2;

        hash ^= k;
        hash = ROT32(hash, r2) * m + n;
    }

    const uint8_t *tail = (const uint8_t *) (key + nblocks * 4);
    uint32_t k1 = 0;

    switch (len & 3) {
    case 3:
        k1 ^= tail[2] << 16;
    case 2:
        k1 ^= tail[1] << 8;
    case 1:
        k1 ^= tail[0];

        k1 *= c1;
        k1 = ROT32(k1, r1);
        k1 *= c2;
        hash ^= k1;
    }

    hash ^= len;
    hash ^= (hash >> 16);
    hash *= 0x85ebca6b;
    hash ^= (hash >> 13);
    hash *= 0xc2b2ae35;
    hash ^= (hash >> 16);

    return hash;
}

uint32_t rebx_hash(const char* str){
    const int reb_seed = 1983;
    return rebx_murmur3_32(str,(uint32_t)strlen(str),reb_seed);
}

/*********************************************************************************
 General particle parameter getter
 ********************************************************************************/

void* rebx_get_param(const struct reb_particle* p, uint32_t param_type){
    struct rebx_param* current = p->ap;
    while(current != NULL){
        if(current->param_type == param_type){
            return current->paramPtr;
        }
        current = current->next;
    }
    return NULL;
}

/*********************************************************************************
 Getters and Setters for particle parameters (need new set for each variable type)
 ********************************************************************************/

void rebx_set_param_double(struct reb_particle* p, const char* param_name, double value){
    uint32_t h = rebx_hash(param_name);
    rebx_set_param_double_hash(p, h, value);
}

void rebx_set_param_double_hash(struct reb_particle* p, uint32_t h, double value){
    double* ptr = rebx_get_param(p, h);
    if(ptr == NULL){
        rebx_add_param_double(p, h, value); // add a new parameter with value if it doesn't exist already
    }else{
        *ptr = value;                       // update existing value
    }
}

void rebx_add_param_double(struct reb_particle* p, uint32_t param_type, double value){
    struct rebx_param* newparam = malloc(sizeof(*newparam));
    newparam->paramPtr = malloc(sizeof(double));
    *(double*) newparam->paramPtr = value;
    newparam->param_type = param_type;

    newparam->next = p->ap;
    p->ap = newparam;

    rebx_add_param_to_be_freed(p->sim->extras, newparam);
}   

double rebx_get_param_double(struct reb_particle* p, const char* param_name){
    uint32_t h = rebx_hash(param_name);
    return rebx_get_param_double_hash(p, h);
}

double rebx_get_param_double_hash(struct reb_particle* p, uint32_t h){
    double* ptr = rebx_get_param(p, h);
    if(ptr == NULL){
        return NAN;
    }else{
        return *ptr;
    }
}

/****************************************
Custom Effect Adders
*****************************************/

void rebx_add_custom_post_timestep_modification(struct rebx_extras* rebx, void (*custom_post_timestep_modification)(struct reb_simulation* const sim, struct rebx_effect* custom_effect), void* custom_params){
    rebx_add_post_timestep_modification(rebx, custom_params, "custom_post_timestep_modification", custom_post_timestep_modification);
}

void rebx_add_custom_force(struct rebx_extras* rebx, void (*custom_force)(struct reb_simulation* const sim, struct rebx_effect* custom_effect), int force_is_velocity_dependent, void* custom_params){
    rebx_add_force(rebx, custom_params, "custom_force", custom_force, force_is_velocity_dependent);
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
