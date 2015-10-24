/**
 * @file 	radiation_forces.c
 * @brief 	Add radiation forces
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
#include <math.h>
#include <stdlib.h>
#include "reboundx.h"

void rebx_radiation_forces(struct reb_simulation* const sim){
	struct reb_particle* particles = sim->particles;
	const struct reb_particle star = particles[0];
	const int N = sim->N;
	struct rebx_extras* rebx = sim->extras;
	const double c = rebx->radiation_forces.c;
	const double L = rebx->radiation_forces.L;
#pragma omp parallel for
	for (int i=1;i<N;i++){
		const struct reb_particle p = particles[i];
		double* Q_pr = rebx_search_param(&p, Q_PR);
		if(Q_pr == NULL) continue; // only particles with Q_pr set feel radiation forces
		const double dx = p.x - star.x; 
		const double dy = p.y - star.y;
		const double dz = p.z - star.z;
		const double dr = sqrt(dx*dx + dy*dy + dz*dz); // distance to star
		
		const double dvx = p.vx - star.vx;
		const double dvy = p.vy - star.vy;
		const double dvz = p.vz - star.vz;
		const double rdot = (dx*dvx + dy*dvy + dz*dvz)/dr; // radial velocity
		const double F_rad = L*p.r*p.r*(*Q_pr)/(4.*c*dr*dr*p.m);

		// Equation (5) of Burns, Lamy & Soter (1979)

		particles[i].ax += F_rad*((1.-rdot/c)*dx/dr - dvx/c);
		particles[i].ay += F_rad*((1.-rdot/c)*dy/dr - dvy/c);
		particles[i].az += F_rad*((1.-rdot/c)*dz/dr - dvz/c);
	}
}
