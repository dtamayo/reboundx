/**
 * @file 	tides.c
 * @brief 	Add tidal forces
 * @author 	Pablo Galindo Salgado, Francisco Pozuelos, Daniel Tamayo  <tamayo.daniel@gmail.com>
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Pablo Galindo Salgado, Francisco Pozuelos, Daniel Tamayo, Hanno Rein
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
#include <math.h>
#include <stdlib.h>
#include "reboundx.h"

void rebx_tides(struct reb_simulation* const sim){
	struct reb_particle* particles = sim->particles;
	const struct rebx_extras* rebx = sim->extras;
	const int _N_tides_active = rebx->tides.N_tides_active;
	const int _N_real = sim->N - sim->N_var; // number of particles in the simulation (took out variational particles)

	for(int i=0; i<_N_tides_active; i++){
		const double k2_1 = rebx_get_tidal_k2(&particles[i]);
		const double tau_1 = rebx_get_tidal_tau(&particles[i]);
		const double M_1 = particles[i].m;
		const double R_1 = particles[i].r;
		const struct rebx_vec3d Omega_1 = rebx_get_rot_Omega(&particles[i]);
		const struct rebx_mom_of_inertia I_1 = rebx_get_mom_of_inertia(&particles[i]);

		for(int j=i+1; j<_N_real; j++){
			const double k2_2 = rebx_get_tidal_k2(&particles[j]);
			const double tau_2 = rebx_get_tidal_tau(&particles[j]);
			const double M_2 = particles[j].m;
			const double R_2 = particles[j].r;
			const struct rebx_vec3d Omega_2 = rebx_get_rot_Omega(&particles[j]);
			const struct rebx_mom_of_inertia I_2 = rebx_get_mom_of_inertia(&particles[j]);
		
			// Do calculation here
		}
	}
}
