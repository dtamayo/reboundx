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
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

void rebx_modify_orbits_direct(struct reb_simulation* const sim){
	struct rebx_extras* rebx = sim->extras;
	struct rebx_params_modify_orbits modparams = rebx->modify_orbits_direct;
	
	struct reb_particle com = {0};
	switch(modparams.coordinates){
	case JACOBI:
		rebxtools_get_com(sim, sim->N-1, &com); // We start with outermost particle, so get COM for the first N-1 particles
		break;
	case BARYCENTRIC:
		rebxtools_get_com(sim, sim->N, &com); // COM of whole system
		break;
	case HELIOCENTRIC:
		com = sim->particles[0];
		break;
	default:
		fprintf(stderr, "Coordinate system not set in modify_orbits_direct?! \n");
		exit(1);
	}

	for(int i=sim->N-1;i>0;--i){
		struct reb_particle* p = &(sim->particles[i]);
		struct reb_orbit o = rebxtools_particle_to_orbit(sim->G, p, &com);
		const double dt = sim->dt_last_done;
		const double tau_a = rebx_get_tau_a(p);
		const double tau_e = rebx_get_tau_e(p);
		const double tau_inc = rebx_get_tau_inc(p);
		const double tau_omega = rebx_get_tau_omega(p);
		const double tau_Omega = rebx_get_tau_Omega(p);

		o.a = o.a*exp(dt/tau_a);
		o.e = o.e*exp(dt/tau_e);
		o.inc = o.inc*exp(dt/tau_inc);
		o.omega += 2*M_PI*dt/tau_omega;
		o.Omega += 2*M_PI*dt/tau_Omega;

		o.a += 2.*o.a*o.e*o.e*modparams.p*dt/tau_e; // Coupling term between e and a

		rebxtools_orbit2p(sim->G, p, &com, &o);
		if(modparams.coordinates == JACOBI){
			rebxtools_update_com_without_particle(&com, p);
		}
	}
	//rebxtools_move_to_com(sim);
}

