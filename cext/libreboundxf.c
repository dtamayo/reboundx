#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "libreboundxf.h"
#include "xftools.h"

//disk parameters for precession
/*double gam;
double Rc;
double diskmass;
double alpha_over_rGM0;
double podot; // pericenter precession at r = Rc
*/
// pointers for damping timescales

void rebxf_add_modify_elements_forces(struct reb_simulation* sim);
void rebxf_add_modify_elements_direct(struct reb_simulation* sim);
void rebxf_add_gr(struct reb_simulation* sim);

struct rebxf_params* rebxf_init(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*) malloc(sizeof(struct rebxf_params));
	xf->tau_a = NULL;
	xf->tau_e = NULL;
	xf->tau_inc = NULL;
	xf->tau_pomega = NULL;

	xf->e_damping_p = 0.;

	xf->allocatedN = 0; 

	xf->forces = NULL;
	xf->Nforces = 0;
	xf->ptm = NULL;
	xf->Nptm = 0;

	sim->xf_params = (struct rebxf_params*) xf;
	return (struct rebxf_params*)sim->xf_params;
}

void rebxf_add(struct reb_simulation* sim, enum REBXFS perturbation){
	switch(perturbation){
		case REBXF_MODIFY_ELEMENTS_FORCES:
			rebxf_add_modify_elements_forces(sim);
			break;
		case REBXF_MODIFY_ELEMENTS_DIRECT:
			rebxf_add_modify_elements_direct(sim);
			break;
		case REBXF_GR:
			rebxf_add_gr(sim);
			break;
		default:
			fprintf(stderr, "Perturbation passed to rebxf_add not recognized.\n");
			exit(1);
	}
}

void rebxf_add_element_timescales(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	
	xf->allocatedN = sim->N;
	xf->tau_a = calloc(xf->allocatedN, sizeof(double));
	xf->tau_e = calloc(xf->allocatedN, sizeof(double));
	xf->tau_inc = calloc(xf->allocatedN, sizeof(double));
	xf->tau_pomega = calloc(xf->allocatedN, sizeof(double));
}
void test(void){
	printf("Hi\n");
}
void rebxf_add_modify_elements_forces(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	rebxf_add_element_timescales(sim);
	xf->Nforces++;
	xf->forces = malloc(sizeof(void*));
	//xf->forces = realloc(xf->forces, sizeof(void*)*xf->Nforces);
	xf->forces[xf->Nforces-1] = test;
	xf->forces[xf->Nforces-1];
}

void rebxf_add_modify_elements_direct(struct reb_simulation* sim){
	rebxf_add_element_timescales(sim);
}

void rebxf_add_gr(struct reb_simulation* sim){
	printf("gr\n");
}

void rebxf_free_xfparams(struct rebxf_params* const xf){
	free(xf->tau_a);
	free(xf->tau_e);
	free(xf->tau_inc);
	free(xf->tau_pomega);
}

void rebxf_check_N(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	if(xf->allocatedN < sim->N) {
		xf->allocatedN = sim->N;
		xf->tau_a = realloc(xf->tau_a, sizeof(double)*xf->allocatedN);
		xf->tau_e = realloc(xf->tau_e, sizeof(double)*xf->allocatedN);
		xf->tau_inc = realloc(xf->tau_inc, sizeof(double)*xf->allocatedN);
		xf->tau_pomega = realloc(xf->tau_pomega, sizeof(double)*xf->allocatedN);
	}
}

double* rebxf_get_tau_a(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	return xf->tau_a;	
}

void rebxf_set_tau_a(struct reb_simulation* sim, double* tau_a){
	rebxf_check_N(sim);
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	memcpy(xf->tau_a, tau_a, sizeof(double)*sim->N);	
}

double* rebxf_get_tau_e(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	return xf->tau_e;	
}

void rebxf_set_tau_e(struct reb_simulation* sim, double* tau_e){
	rebxf_check_N(sim);
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	memcpy(xf->tau_e, tau_e, sizeof(double)*sim->N);	
}

double* rebxf_get_tau_inc(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	return xf->tau_inc;	
}

void rebxf_set_tau_inc(struct reb_simulation* sim, double* tau_inc){
	rebxf_check_N(sim);
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	memcpy(xf->tau_inc, tau_inc, sizeof(double)*sim->N);	
}

double* rebxf_get_tau_pomega(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	return xf->tau_pomega;	
}

void rebxf_set_tau_pomega(struct reb_simulation* sim, double* tau_pomega){
	rebxf_check_N(sim);
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	memcpy(xf->tau_pomega, tau_pomega, sizeof(double)*sim->N);	
}
