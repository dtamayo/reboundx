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
#include "reboundx.h"

//disk parameters for precession
/*double gam;
double Rc;
double diskmass;
double alpha_over_rGM0;
double podot; // pericenter precession at r = Rc
*/
// pointers for damping timescales

struct rebx_extras* rebx_init(struct reb_simulation* sim){
	struct rebx_extras* rebx = malloc(sizeof(struct rebx_extras));
	rebx_initialize(sim, rebx);
	return rebx;
}

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx){
	sim->extras = rebx;
	rebx->sim = sim;

	rebx->forces = NULL;
	rebx->Nforces = 0;
	rebx->ptm = NULL;
	rebx->Nptm = 0;
	
	rebx->modify_orbits_forces->e_damping_p = 0.;
	rebx->modify_orbits_direct->e_damping_p = 0.;

	rebx->gr->c = C_DEFAULT; // speed of light in default units of AU/(yr/2pi)
}

/*void rebx_add(struct reb_simulation* sim, enum REBX_EXTRAS extra){
	switch(extra){
		case REBX_MODIFY_ORBITS_FORCES:
			rebx_add_modify_orbits_forces(sim);
			break;
		case REBX_MODIFY_ORBITS_DIRECT:
			rebx_add_modify_orbits_direct(sim);
			break;
		case REBX_GR:
			rebx_add_gr(sim, C_DEFAULT);
			break;
		default:
			fprintf(stderr, "Perturbation passed to rebx_add not recognized.\n");
			exit(1);
	}
}*/

void rebx_forces(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	for(int i=sim->N-1; i>=0; i--){
		struct rebx_p_param* current = sim->particles[i].ap;
		while(current != NULL){

		}
	}
	/*for(int i=0;i<rebx->Nforces;i++){
		rebx->forces[i](sim);
	}*/
}

void rebx_ptm(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	for(int i=0;i<rebx->Nptm;i++){
		rebx->ptm[i](sim);
	}
}

void rebx_add_modify_orbits_forces(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	
	rebx->Nforces++;
	if(!rebx->forces){
		rebx->forces = malloc(sizeof(xptr));
	}
	else{
		rebx->forces = realloc(rebx->forces, sizeof(xptr)*rebx->Nforces);
	}
	rebx->forces[rebx->Nforces-1] = rebx_modify_orbits_forces;
	sim->additional_forces = rebx_forces;
	sim->force_is_velocity_dependent = 1;
}

void rebx_add_modify_orbits_direct(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	
	rebx->Nptm++;
	if(!rebx->ptm){
		rebx->ptm = malloc(sizeof(xptr));
	}
	else{
		rebx->ptm = realloc(rebx->ptm, sizeof(xptr)*rebx->Nptm);
	}
	rebx->ptm[rebx->Nptm-1] = rebx_modify_orbits_direct;
	sim->post_timestep_modifications = rebx_ptm;
	//sim->pre_timestep_modifications = rebx_ptm;
}

void rebx_add_gr(struct reb_simulation* sim, double c){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	rebx->gr->c = c;
	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(xptr)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_gr;
	sim->additional_forces = rebx_forces;
	sim->force_is_velocity_dependent = 1;
}

void rebx_add_gr_single_mass(struct reb_simulation* sim, double c){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	rebx->gr->c = c;
	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(xptr)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_gr_single_mass;
	sim->additional_forces = rebx_forces;
	sim->force_is_velocity_dependent = 1;
}

void rebx_add_gr_potential(struct reb_simulation* sim, double c){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	rebx->gr->c = c;
	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(xptr)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_gr_potential;
	sim->additional_forces = rebx_forces;
}

/* Garbage collection */
void rebx_free(struct rebx_extras* const rebx){
	rebx_free_pointers(rebx);
	free(rebx);
}

void rebx_free_pointers(struct rebx_extras* const rebx){
	if(rebx->forces){
		free(rebx->forces);
	}
	if(rebx->ptm){
		free(rebx->ptm);
	}
}

/*Getter/setters*/

/*void rebx_check_N(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	if(rebx->allocatedN < sim->N) {
		rebx->allocatedN = sim->N;
		rebx->tau_a = realloc(rebx->tau_a, sizeof(double)*rebx->allocatedN);
		rebx->tau_e = realloc(rebx->tau_e, sizeof(double)*rebx->allocatedN);
		rebx->tau_inc = realloc(rebx->tau_inc, sizeof(double)*rebx->allocatedN);
		rebx->tau_omega = realloc(rebx->tau_omega, sizeof(double)*rebx->allocatedN);
	}
}*/

/*double* rebx_get_tau_a(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	return rebx->tau_a;	
}

void rebx_set_tau_a(struct reb_simulation* sim, double* tau_a){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	memcpy(rebx->tau_a, tau_a, sizeof(double)*sim->N);	
}

double* rebx_get_tau_e(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	return rebx->tau_e;	
}

void rebx_set_tau_e(struct reb_simulation* sim, double* tau_e){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	memcpy(rebx->tau_e, tau_e, sizeof(double)*sim->N);	
}

double* rebx_get_tau_inc(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	return rebx->tau_inc;	
}

void rebx_set_tau_inc(struct reb_simulation* sim, double* tau_inc){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	memcpy(rebx->tau_inc, tau_inc, sizeof(double)*sim->N);	
}

double* rebx_get_tau_omega(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	return rebx->tau_omega;	
}

void rebx_set_tau_omega(struct reb_simulation* sim, double* tau_omega){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	memcpy(rebx->tau_omega, tau_omega, sizeof(double)*sim->N);	
}*/
