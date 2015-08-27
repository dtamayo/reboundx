#include "elements_direct.h"
#include "reboundxf.h"
#include "xftools.h"
#include <stdio.h>
#include <stdlib.h>

void rebxf_elements_direct(struct reb_simulation* const sim){
	//rebxf_check_N(sim);
	struct rebxf_params* xf = sim->xf_params;
	struct rebxf_param_elements_direct* xfparams = &xf->elements_direct;
	struct reb_particle com = sim->particles[0];
	for(int i=1;i<sim->N;i++){
		struct reb_particle *p = &(sim->particles[i]);
		int* err = malloc(sizeof(int)); // dummy
		struct reb_orbit o = xftools_particle_to_orbit(sim->G, sim->particles[i], com, err);
	    double da = 0.;
		double de = 0.;
		double dom = 0.;	
		if (xfparams->tau_a[i] != 0.){
			da += -o.a*sim->dt/xfparams->tau_a[i]; 
		}
	
		if (xfparams->tau_e[i] != 0.){
			de += -o.e*sim->dt/xfparams->tau_e[i];
			da += -2.*o.a*o.e*o.e*xfparams->e_damping_p*sim->dt/xfparams->tau_e[i];
		}

		if (xfparams->tau_omega[i] != 0.){
			dom += 2*M_PI*sim->dt/xfparams->tau_omega[i];
		}

		o.a += da;
		o.e += de;
		o.omega += dom;

		xftools_orbit2p(sim->G, &sim->particles[i], &com, o); 
	}
	xftools_move_to_com(sim);
}

