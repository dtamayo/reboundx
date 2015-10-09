/**
 * @file 	modify_orbits_direct.c
 * @brief 	Update orbital elements at the end of each timestep
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "modify_orbits_direct.h"
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

void rebx_modify_orbits_direct(struct reb_simulation* const sim){
	//rebx_check_N(sim);
	struct rebx_params_modify_orbits_direct* rebxparams = &((struct rebx_extras*)(sim->extras))->modify_orbits_direct;
	struct reb_particle com = rebxtools_get_com(sim);
	
	for(int i=sim->N-1;i>0;--i){
		struct reb_particle *p = &(sim->particles[i]);
		com = rebxtools_get_com_without_particle(com, *p);
		int* err = malloc(sizeof(int)); // dummy
		struct reb_orbit o = rebxtools_particle_to_orbit_err(sim->G, sim->particles[i], com, err);
	    double da = 0.;
		double de = 0.;
		double di = 0.;
		double dom = 0.;	
		double dt2 = sim->dt_last_done;
		if (rebxparams->tau_a[i] != 0.){
			da += o.a*sim->dt/rebxparams->tau_a[i]; 
		}
	
		if (rebxparams->tau_e[i] != 0.){
			o.e = o.e*exp(dt2/rebxparams->tau_e[i]);
			//o.e += o.e*(exp(sim->dt/rebxparams->tau_e[i])-1.);
			//de += o.e*sim->dt/rebxparams->tau_e[i];
			da += 2.*o.a*o.e*o.e*rebxparams->e_damping_p*sim->dt/rebxparams->tau_e[i];
		}

		if (rebxparams->tau_inc[i] != 0.){
			o.inc = o.inc*exp(dt2/rebxparams->tau_inc[i]);
			//di += o.inc*sim->dt/rebxparams->tau_inc[i];
		}

		if (rebxparams->tau_omega[i] != 0.){
			dom += 2*M_PI*sim->dt/rebxparams->tau_omega[i];
		}

		o.a += da;
		o.e += de;
		o.inc += di;
		o.omega += dom;

		rebxtools_orbit2p(sim->G, &sim->particles[i], &com, o); 
	}
	//rebxtools_move_to_com(sim);
}

