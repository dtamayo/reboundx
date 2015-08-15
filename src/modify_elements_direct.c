#include "modify_elements_direct.h"
#include <stdio.h>
void rebxf_modify_elements_direct(struct reb_simulation* const sim){
	printf("rebxf_modify_elements_direct\n");
	/*rebxf_check_N(sim);
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	struct reb_particle com = xftools_get_com(sim);
	for(int i=1;i<sim->N;i++){
		struct reb_particle *p = &(sim->particles[i]);
		struct reb_orbit o = xftools_p2orbit(sim->G, sim->particles[i], com);
	    double da = 0.;
		double de = 0.;
		double dpo = 0.;	
		if (xf->tau_a[i] != 0.){
			da += -o.a*sim->dt/xf->tau_a[i]; 
		}
		
		if (xf->tau_e[i] != 0.){
			de += -o.e*sim->dt/xf->tau_e[i];
			da += -2.*o.a*o.e*o.e*xf->e_damping_p*sim->dt/xf->tau_e[i];
		}

		if (xf->tau_pomega[i] != 0.){
			dpo += 2*M_PI*sim->dt/xf->tau_pomega[i]*(1.+sin(o.omega));
		}

		o.a += da;
		o.e += de;
		o.omega += dpo;

		xftools_orbit2p(&sim->particles[i], sim->G, &com, o); 
	}
	xftools_move_to_com(sim->particles, sim->N);*/
}

