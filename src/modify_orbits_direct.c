#include "modify_orbits_direct.h"
#include "reboundx.h"
#include "rebxtools.h"
#include <stdio.h>
#include <stdlib.h>

void rebx_modify_orbits_direct(struct reb_simulation* const sim){
	//rebx_check_N(sim);
	struct rebx_params_modify_orbits_direct* rebxparams = &((struct rebx_extras*)(sim->extras))->modify_orbits_direct;
	struct reb_particle com = sim->particles[0];
	for(int i=1;i<sim->N;i++){
		struct reb_particle *p = &(sim->particles[i]);
		int* err = malloc(sizeof(int)); // dummy
		struct reb_orbit o = rebxtools_particle_to_orbit(sim->G, sim->particles[i], com, err);
	    double da = 0.;
		double de = 0.;
		double dom = 0.;	
		if (rebxparams->tau_a[i] != 0.){
			da += -o.a*sim->dt/rebxparams->tau_a[i]; 
		}
	
		if (rebxparams->tau_e[i] != 0.){
			de += -o.e*sim->dt/rebxparams->tau_e[i];
			da += -2.*o.a*o.e*o.e*rebxparams->e_damping_p*sim->dt/rebxparams->tau_e[i];
		}

		if (rebxparams->tau_omega[i] != 0.){
			dom += 2*M_PI*sim->dt/rebxparams->tau_omega[i];
		}

		o.a += da;
		o.e += de;
		o.omega += dom;

		rebxtools_orbit2p(sim->G, &sim->particles[i], &com, o); 
	}
	rebxtools_move_to_com(sim);
}

