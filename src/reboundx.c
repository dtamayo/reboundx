/**
 * @file 	reboundx.c
 * @brief 	Initialization and getter/setter routines.
 * @author 	Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section 	LICENSE
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

//disk parameters for precession
/*double gam;
double Rc;
double diskmass;
double alpha_over_rGM0;
double podot; // pericenter precession at r = Rc
*/
// pointers for damping timescales

void rebx_set_double(struct reb_simulation* sim, int p_index, enum REBX_PARAMS param, double value){
	struct reb_particle* p = &(sim->particles[p_index]);
	rebx_add_double_param(sim, &(p->ap), param, value);
}

void rebx_add_particle(struct rebx_extras* rebx, struct rebx_p_param* p_param){
	while(rebx->allocatedN <= rebx->N){
		rebx->allocatedN += 128;
		rebx->particles = realloc(rebx->particles, sizeof(struct rebx_p_param*)*rebx->allocatedN);
	}

	p_param->rebx_index = rebx->N;	// uniquely identifies this linked list in rebx->particles
	(rebx->N)++;
}

void rebx_update_particles(struct rebx_extras* rebx, struct rebx_p_param* p_param){
	if(p_param->next == NULL){	// This was the first parameter added to particle, so need to add to rebx->particles
		rebx_add_particle(rebx, p_param);
	}
	else{
		p_param->rebx_index = p_param->next->rebx_index;
	}
	rebx->particles[p_param->rebx_index] = p_param;
}

void rebx_add_double_param(struct reb_simulation* sim, void** _p_paramsRef, enum REBX_PARAMS param, double value){
	struct rebx_p_param** p_paramsRef = (struct rebx_p_param**)_p_paramsRef; // avoid void** to p_param** warning
	struct rebx_p_param* newparam = malloc(sizeof(struct rebx_p_param));
	
	newparam->valPtr = malloc(sizeof(double));
	*(double*)newparam->valPtr = value; // valPtr is a void* so need to cast before dereference
	newparam->param = param;
	newparam->next = *p_paramsRef;
	*p_paramsRef = newparam;

	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	rebx_update_particles(rebx, *p_paramsRef);
}

void rebx_free_p_params(struct rebx_p_param* apPtr){
	if(apPtr != NULL){
		rebx_free_p_params(apPtr->next);
		free(apPtr->valPtr);
		free(apPtr->next);
	}
}

void rebx_free_particles(struct rebx_extras* rebx){
	for(int i=0; i<rebx->N; i++)
	{
		rebx_free_p_params(rebx->particles[i]);
		free(rebx->particles[i]);
	}
}

void rebx_free(struct rebx_extras* rebx){
	rebx_free_particles(rebx);
	free(rebx->particles);
	free(rebx);
}

double rebx_get_double(struct reb_particle p, enum REBX_PARAMS param){
	struct rebx_p_param* current;
	for(current = p.ap; current != NULL; current=current->next){
		if(current->param == param){
			return *(double*)current->valPtr;
		}
	}
	fprintf(stderr, "REBOUNDx parameter not found in call to rebx_get_double.  Returning 0.\n");
	return 0.;
}

double rebx_get_double_param(struct rebx_p_param* p_param){
	return *(double*)p_param->valPtr;
}

void rebx_add_int_param(struct rebx_p_param** p_paramsRef, enum REBX_PARAMS param, int value){
	struct rebx_p_param* newparam = malloc(sizeof(struct rebx_p_param));

	newparam->valPtr = malloc(sizeof(int));
	*(int*)newparam->valPtr = value;
	newparam->param = param;
	newparam->next = *p_paramsRef;
	*p_paramsRef = newparam;
}

int rebx_get_int_param(struct rebx_p_param* p_param){
	return *(int*)p_param->valPtr;
}

struct rebx_extras* rebx_init(struct reb_simulation* sim){
	struct rebx_extras* rebx = malloc(sizeof(struct rebx_extras));
	rebx_initialize(sim, rebx);
	return rebx;
}

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx){
	sim->extras = rebx;
	rebx->sim = sim;

	rebx->N = 0;
	rebx->allocatedN = 0;
	rebx->particles = NULL;

	rebx->modify_orbits_forces.e_damping_p = 0.;
	rebx->modify_orbits_forces.coordinates = JACOBI;

	rebx->modify_orbits_direct.e_damping_p = 0.;
	rebx->modify_orbits_direct.coordinates = JACOBI;

	rebx->gr.c = C_DEFAULT; // speed of light in default units of AU/(yr/2pi)

	rebx->ptm = NULL;
	rebx->forces = NULL;
	rebx->Nptm = 0;
	rebx->Nforces = 0;
}

void rebx_forces(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	for(int i=0; i<rebx->Nforces; i++){
		rebx->forces[i](sim);
	}
}

void rebx_ptm(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	for(int i=0; i<rebx->Nforces; i++){
		rebx->ptm[i](sim);
	}
}

void rebx_add_gr(struct reb_simulation* sim, double c){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	sim->additional_forces = rebx_forces;
	sim->force_is_velocity_dependent = 1;

	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_gr;
	rebx->gr.c = c;
}

void rebx_add_gr_single_mass(struct reb_simulation* sim, double c){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	sim->additional_forces = rebx_forces;
	sim->force_is_velocity_dependent = 1;

	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_gr_single_mass;
	rebx->gr.c = c;
}

void rebx_add_gr_potential(struct reb_simulation* sim, double c){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	sim->additional_forces = rebx_forces;

	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_gr_potential;
	rebx->gr.c = c;
}

void rebx_add_modify_orbits_forces(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	sim->additional_forces = rebx_forces;
	sim->force_is_velocity_dependent = 1;

	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_modify_orbits_forces;
}

void rebx_add_modify_orbits_direct(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	sim->post_timestep_modifications = rebx_ptm;

	rebx->Nptm++;
	rebx->ptm = realloc(rebx->ptm, sizeof(*rebx->ptm)*rebx->Nptm);
	rebx->ptm[rebx->Nptm-1] = rebx_modify_orbits_direct;
}
