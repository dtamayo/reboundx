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
		printf("t=%e\tdt=%e\tptm\n", sim->t, sim->dt);
		rebx->ptm[i](sim);
	}
}

/* Initialization routines. */

struct rebx_extras* rebx_init(struct reb_simulation* sim){
	struct rebx_extras* rebx = malloc(sizeof(*rebx));
	rebx_initialize(sim, rebx);
	return rebx;
}

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx){
	sim->extras = rebx;
	rebx->sim = sim;

	rebx->params_to_be_freed = NULL;
	rebx->ptm = NULL;
	rebx->forces = NULL;

	rebx->modify_orbits_forces = NULL;
	rebx->modify_orbits_direct = NULL;
	rebx->gr = NULL;

	rebx->Nptm = 0;
	rebx->Nforces = 0;
}

/* Garbage collection routines. */

void rebx_free_params(struct rebx_extras* rebx){
	struct rebx_param_to_be_freed* current = rebx->params_to_be_freed;
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

/* Internal utility functions. */

void rebx_add_param_to_be_freed(struct rebx_param_to_be_freed** ptbfRef, struct rebx_param* param){
	struct rebx_param_to_be_freed* newparam = malloc(sizeof(*newparam));
	newparam->param = param;

	newparam->next = *ptbfRef;
	*ptbfRef = newparam;
}

void* rebx_search_param(struct rebx_param* current, enum REBX_PARAMS param){
	while(current != NULL){
		if(current->param == param){
			return current;
		}
	}
	return NULL;
}

/* Internal parameter adders (need a different one for each REBX_PARAM type). */

void rebx_add_param_orb_tau(struct reb_simulation* sim, void** paramsRef){
	struct rebx_param* newparam = malloc(sizeof(*newparam));

	newparam->valPtr = malloc(sizeof(struct rebx_orb_tau));
	struct rebx_orb_tau orb_tau = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY};
	*(struct rebx_orb_tau*) newparam->valPtr = orb_tau;
	newparam->param = (enum REBX_PARAMS)ORB_TAU;
	newparam->next = *paramsRef;
	*paramsRef = newparam;

	struct rebx_extras* rebx = sim->extras;
	rebx_add_param_to_be_freed(&rebx->params_to_be_freed, *paramsRef);
}

/* User-called getters and setters for each parameter*/

void rebx_set_tau_a(struct reb_simulation* sim, int p_index, double value){
	struct reb_particle* p = &(sim->particles[p_index]);
	struct rebx_param* orb_tau = rebx_search_param(p->ap, ORB_TAU);
	if(orb_tau == NULL){
		rebx_add_param_orb_tau(sim, &p->ap);
		void* valPtr = ((struct rebx_param*) p->ap)->valPtr;
		((struct rebx_orb_tau*) valPtr)->tau_a = value;
	}
	else{
		((struct rebx_orb_tau*) (orb_tau->valPtr))->tau_a = value;
	}
}

void rebx_set_tau_e(struct reb_simulation* sim, int p_index, double value){
	struct reb_particle* p = &(sim->particles[p_index]);
	struct rebx_param* orb_tau = rebx_search_param(p->ap, ORB_TAU);
	if(orb_tau == NULL){
		rebx_add_param_orb_tau(sim, &p->ap);
		void* valPtr = ((struct rebx_param*) p->ap)->valPtr;
		((struct rebx_orb_tau*) valPtr)->tau_e = value;
	}
	else{
		((struct rebx_orb_tau*) (orb_tau->valPtr))->tau_e = value;
	}
}

void rebx_set_tau_inc(struct reb_simulation* sim, int p_index, double value){
	struct reb_particle* p = &(sim->particles[p_index]);
	struct rebx_param* orb_tau = rebx_search_param(p->ap, ORB_TAU);
	if(orb_tau == NULL){
		rebx_add_param_orb_tau(sim, &p->ap);
		void* valPtr = ((struct rebx_param*) p->ap)->valPtr;
		((struct rebx_orb_tau*) valPtr)->tau_inc = value;
	}
	else{
		((struct rebx_orb_tau*) (orb_tau->valPtr))->tau_inc = value;
	}
}

void rebx_set_tau_omega(struct reb_simulation* sim, int p_index, double value){
	struct reb_particle* p = &(sim->particles[p_index]);
	struct rebx_param* orb_tau = rebx_search_param(p->ap, ORB_TAU);
	if(orb_tau == NULL){
		rebx_add_param_orb_tau(sim, &p->ap);
		void* valPtr = ((struct rebx_param*) p->ap)->valPtr;
		((struct rebx_orb_tau*) valPtr)->tau_omega = value;
	}
	else{
		((struct rebx_orb_tau*) (orb_tau->valPtr))->tau_omega = value;
	}
}

void rebx_set_tau_Omega(struct reb_simulation* sim, int p_index, double value){
	struct reb_particle* p = &(sim->particles[p_index]);
	struct rebx_param* orb_tau = rebx_search_param(p->ap, ORB_TAU);
	if(orb_tau == NULL){
		rebx_add_param_orb_tau(sim, &p->ap);
		void* valPtr = ((struct rebx_param*) p->ap)->valPtr;
		((struct rebx_orb_tau*) valPtr)->tau_Omega = value;
	}
	else{
		((struct rebx_orb_tau*) (orb_tau->valPtr))->tau_Omega = value;
	}
}

double rebx_get_tau_a(struct reb_particle p){
	struct rebx_param* orb_tau = rebx_search_param(p.ap, ORB_TAU);
	if(orb_tau == NULL){
		fprintf(stderr, "tau_a wasn't set for particle passed to rebx_get_tau_a\n");
	}
	return ((struct rebx_orb_tau*) (orb_tau->valPtr))->tau_a;
}

double rebx_get_tau_e(struct reb_particle p){
	struct rebx_param* orb_tau = rebx_search_param(p.ap, ORB_TAU);
	if(orb_tau == NULL){
		fprintf(stderr, "tau_e wasn't set for particle passed to rebx_get_tau_e\n");
	}
	return ((struct rebx_orb_tau*) (orb_tau->valPtr))->tau_e;
}

double rebx_get_tau_inc(struct reb_particle p){
	struct rebx_param* orb_tau = rebx_search_param(p.ap, ORB_TAU);
	if(orb_tau == NULL){
		fprintf(stderr, "tau_inc wasn't set for particle passed to rebx_get_tau_inc\n");
	}
	return ((struct rebx_orb_tau*) (orb_tau->valPtr))->tau_inc;
}

double rebx_get_tau_omega(struct reb_particle p){
	struct rebx_param* orb_tau = rebx_search_param(p.ap, ORB_TAU);
	if(orb_tau == NULL){
		fprintf(stderr, "tau_omega wasn't set for particle passed to rebx_get_tau_omega\n");
	}
	return ((struct rebx_orb_tau*) (orb_tau->valPtr))->tau_omega;
}

double rebx_get_tau_Omega(struct reb_particle p){
	struct rebx_param* orb_tau = rebx_search_param(p.ap, ORB_TAU);
	if(orb_tau == NULL){
		fprintf(stderr, "tau_Omega wasn't set for particle passed to rebx_get_tau_Omega\n");
	}
	return ((struct rebx_orb_tau*) (orb_tau->valPtr))->tau_Omega;
}

/* Functions to add effects. */

void rebx_add_gr(struct reb_simulation* sim, double c){
	sim->additional_forces = rebx_forces;
	sim->force_is_velocity_dependent = 1;

	struct rebx_extras* rebx = sim->extras;
	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_gr;
	rebx->gr = malloc(sizeof(*rebx->gr));
	rebx->gr->c = c;
}

void rebx_add_gr_single_mass(struct reb_simulation* sim, double c){
	sim->additional_forces = rebx_forces;
	sim->force_is_velocity_dependent = 1;

	struct rebx_extras* rebx = sim->extras;
	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_gr_single_mass;
	rebx->gr = malloc(sizeof(*rebx->gr));
	rebx->gr->c = c;
}

void rebx_add_gr_potential(struct reb_simulation* sim, double c){
	sim->additional_forces = rebx_forces;

	struct rebx_extras* rebx = sim->extras;
	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_gr_potential;
	rebx->gr = malloc(sizeof(*rebx->gr));
	rebx->gr->c = c;
}

void rebx_add_modify_orbits_forces(struct reb_simulation* sim){
	sim->additional_forces = rebx_forces;
	sim->force_is_velocity_dependent = 1;

	struct rebx_extras* rebx = sim->extras;
	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(*rebx->forces)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_modify_orbits_forces;
	rebx->modify_orbits_forces = malloc(sizeof(*rebx->modify_orbits_forces));
	rebx->modify_orbits_forces->e_damping_p = 0.;
	rebx->modify_orbits_forces->coordinates = JACOBI;
}

void rebx_add_modify_orbits_direct(struct reb_simulation* sim){
	sim->post_timestep_modifications = rebx_ptm;

	struct rebx_extras* rebx = sim->extras;
	rebx->Nptm++;
	rebx->ptm = realloc(rebx->ptm, sizeof(*rebx->ptm)*rebx->Nptm);
	rebx->ptm[rebx->Nptm-1] = rebx_modify_orbits_direct;
	rebx->modify_orbits_direct = malloc(sizeof(*rebx->modify_orbits_direct));
	rebx->modify_orbits_direct->e_damping_p = 0.;
	rebx->modify_orbits_direct->coordinates = JACOBI;
}
