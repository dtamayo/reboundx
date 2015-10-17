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

void* rebx_search_param(struct rebx_param* current, enum REBX_PARAMS param){
	while(current != NULL){
		if(current->param == param){
			return current;
		}
	}
	return NULL;
}

void rebx_add_param_orb_tau(struct reb_simulation* sim, struct rebx_param** paramsRef){
	struct rebx_param* newparam = malloc(sizeof(struct rebx_param));

	newparam->valPtr = malloc(sizeof(struct rebx_orb_tau));
	newparam->param = (enum REBX_PARAMS)ORB_TAU;
	newparam->next = *paramsRef;
	*paramsRef = newparam;

	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	rebx_add_param_to_be_freed(rebx, *paramsRef);
}

void rebx_set_tau_a(struct reb_simulation* sim, int p_index, double value){
	struct reb_particle* p = &(sim->particles[p_index]);
	struct rebx_param* ap = (struct rebx_param*)p->ap;
	struct rebx_param* orb_tau = rebx_search_param(ap, ORB_TAU);
	if(orb_tau == NULL){
		rebx_add_param_orb_tau(sim, &ap);
		ap = (struct rebx_param*)p->ap;
		((struct rebx_orb_tau*) ap->valPtr)->tau_a = value;
	}
	else{
		((struct rebx_orb_tau*) (orb_tau->valPtr))->tau_a = value;
	}
}

double rebx_get_tau_a(struct reb_particle p){
	struct rebx_param* orb_tau = rebx_search_param(p.ap, ORB_TAU);
	if(orb_tau == NULL){
		fprintf(stderr, "tau_a wasn't set for particle passed to rebx_get_tau_a\n");
	}
	return ((struct rebx_orb_tau*) (orb_tau->valPtr))->tau_a;
}

void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param){
	struct rebx_param_to_be_freed* newparam = malloc(sizeof(*newparam));
	newparam->param = param;
	newparam->next = *rebx->params_to_be_freed;
	*rebx->params_to_be_freed = newparam;
}

/*void rebx_add_double_param(struct reb_simulation* sim, void** _paramsRef, enum REBX_PARAMS param, double value){
	struct rebx_param** paramsRef = (struct rebx_param**)_paramsRef; // avoid void** to param** warning
	struct rebx_param* newparam = malloc(sizeof(struct rebx_param));
	
	newparam->valPtr = malloc(sizeof(double));
	*(double*)newparam->valPtr = value; // valPtr is a void* so need to cast before dereference
	newparam->param = param;
	newparam->next = *paramsRef;
	*paramsRef = newparam;

	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	rebx_update_particles(rebx, *paramsRef);
}*/

void rebx_free_params(struct rebx_extras* rebx){
	struct rebx_param_to_be_freed* current = *rebx->params_to_be_freed;
	struct rebx_param_to_be_freed* temp_next;
	while(current != NULL)
	{
		temp_next = current->next;
		free(current->param->valPtr);
		free(current->param);
		free(current);
		current = temp_next;
	}
}

void rebx_free(struct rebx_extras* rebx){
	rebx_free_params(rebx);
	free(rebx->modify_orbits_forces);
	free(rebx->modify_orbits_direct);
	free(rebx->gr);
	free(rebx->forces);
	free(rebx->ptm);
	free(rebx);
}

struct rebx_extras* rebx_init(struct reb_simulation* sim){
	struct rebx_extras* rebx = malloc(sizeof(struct rebx_extras));
	rebx_initialize(sim, rebx);
	return rebx;
}

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx){
	sim->extras = rebx;
	rebx->sim = sim;

	rebx->params_to_be_freed = NULL;
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
	rebx->gr->c = c;
}

void rebx_add_gr_single_mass(struct reb_simulation* sim, double c){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	sim->additional_forces = rebx_forces;
	sim->force_is_velocity_dependent = 1;

	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_gr_single_mass;
	rebx->gr->c = c;
}

void rebx_add_gr_potential(struct reb_simulation* sim, double c){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	sim->additional_forces = rebx_forces;

	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_gr_potential;
	rebx->gr->c = c;
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
