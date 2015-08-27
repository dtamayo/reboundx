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
	struct rebx_extras* rebx = (struct rebx_extras*) malloc(sizeof(struct rebx_extras));
	sim->extras = (struct rebx_extras*) rebx;
	rebx->sim = sim;

	rebx->forces = NULL;
	rebx->Nforces = 0;
	rebx->ptm = NULL;
	rebx->Nptm = 0;
	
	rebx->modify_orbits_forces.allocatedN = 0; 
	rebx->modify_orbits_forces.tau_a = NULL;
	rebx->modify_orbits_forces.tau_e = NULL;
	rebx->modify_orbits_forces.tau_inc = NULL;
	rebx->modify_orbits_forces.tau_omega = NULL;
	rebx->modify_orbits_forces.e_damping_p = 0.;
	
	rebx->modify_orbits_direct.allocatedN = 0; 
	rebx->modify_orbits_direct.tau_a = NULL;
	rebx->modify_orbits_direct.tau_e = NULL;
	rebx->modify_orbits_direct.tau_inc = NULL;
	rebx->modify_orbits_direct.tau_omega = NULL;
	rebx->modify_orbits_direct.e_damping_p = 0.;

	rebx->gr.c = 10064.9150404; // speed of light in units of AU/(yr/2pi)
	return (struct rebx_extras*)sim->extras;
}

/*void rebx_add(struct reb_simulation* sim, enum REBx_MODS perturbation){
	switch(perturbation){
		case REBx_MODIFY_modify_orbits_forces:
			rebx_add_modify_modify_orbits_forces(sim);
			break;
		case REBx_MODIFY_modify_orbits_direct:
			rebx_add_modify_modify_orbits_direct(sim);
			break;
		case REBx_GR:
			rebx_add_gr(sim);
			break;
		default:
			fprintf(stderr, "Perturbation passed to rebx_add not recognized.\n");
			exit(1);
	}
}*/

void rebx_forces(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	for(int i=0;i<rebx->Nforces;i++){
		rebx->forces[i](sim);
	}
}

void rebx_ptm(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	for(int i=0;i<rebx->Nptm;i++){
		rebx->ptm[i](sim);
	}
}

void rebx_add_modify_orbits_forces(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	rebx->modify_orbits_forces.allocatedN = sim->N;
	rebx->modify_orbits_forces.tau_a = calloc(rebx->modify_orbits_forces.allocatedN, sizeof(double));
	rebx->modify_orbits_forces.tau_e = calloc(rebx->modify_orbits_forces.allocatedN, sizeof(double));
	rebx->modify_orbits_forces.tau_inc = calloc(rebx->modify_orbits_forces.allocatedN, sizeof(double));
	rebx->modify_orbits_forces.tau_omega = calloc(rebx->modify_orbits_forces.allocatedN, sizeof(double));
	
	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(xptr)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_modify_orbits_forces;
	sim->additional_forces = rebx_forces;
	sim->force_is_velocity_dependent = 1;
}

void rebx_add_modify_orbits_direct(struct reb_simulation* sim){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	rebx->modify_orbits_direct.allocatedN = sim->N;
	rebx->modify_orbits_direct.tau_a = calloc(rebx->modify_orbits_direct.allocatedN, sizeof(double));
	rebx->modify_orbits_direct.tau_e = calloc(rebx->modify_orbits_direct.allocatedN, sizeof(double));
	rebx->modify_orbits_direct.tau_inc = calloc(rebx->modify_orbits_direct.allocatedN, sizeof(double));
	rebx->modify_orbits_direct.tau_omega = calloc(rebx->modify_orbits_direct.allocatedN, sizeof(double));
	
	rebx->Nptm++;
	rebx->ptm = realloc(rebx->ptm, sizeof(xptr)*rebx->Nptm);
	rebx->ptm[rebx->Nptm-1] = rebx_modify_orbits_direct;
	sim->post_timestep_modifications = rebx_ptm;
}

void rebx_add_gr(struct reb_simulation* sim, double c){
	struct rebx_extras* rebx = (struct rebx_extras*)sim->extras;
	rebx->gr.c = c;
	rebx->Nforces++;
	rebx->forces = realloc(rebx->forces, sizeof(xptr)*rebx->Nforces);
	rebx->forces[rebx->Nforces-1] = rebx_gr;
	sim->additional_forces = rebx_forces;
	sim->force_is_velocity_dependent = 1;
}

/* Garbage collection */
void rebx_free_xparams(struct rebx_extras* const rebx){
	free(rebx->modify_orbits_forces.tau_a);
	free(rebx->modify_orbits_forces.tau_e);
	free(rebx->modify_orbits_forces.tau_inc);
	free(rebx->modify_orbits_forces.tau_omega);

	free(rebx->modify_orbits_direct.tau_a);
	free(rebx->modify_orbits_direct.tau_e);
	free(rebx->modify_orbits_direct.tau_inc);
	free(rebx->modify_orbits_direct.tau_omega);
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
