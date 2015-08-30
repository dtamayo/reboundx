/**
 * @file 	gr.c
 * @brief 	Post-newtonian general relativity corrections
 * @author 	Pengshuai (Sam) Shi <tamayo.daniel@gmail.com>
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Pengshuai (Sam) Shi, Dan Tamayo, Hanno Rein
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
#include "rebound.h"
#include "reboundx.h"
#include "gr.h"

void rebx_gr(struct reb_simulation* const sim){
	struct rebx_params_gr* rebxparams = &((struct rebx_extras*)(sim->extras))->gr;
	const double C = rebxparams->c;
	const int _N_real = sim->N - sim->N_var;
	const double G = sim->G;
	struct reb_particle* const particles = sim->particles;
	
	const struct reb_particle sun = particles[0];
	for (int i=1; i<_N_real; i++){
		const double dx = particles[i].x - sun.x;
		const double dy = particles[i].y - sun.y;
		const double dz = particles[i].z - sun.z;
		const double r2 = dx*dx + dy*dy + dz*dz;
		const double _r = sqrt(r2);
		const double dvx = particles[i].vx - sun.vx;
		const double dvy = particles[i].vy - sun.vy;
		const double dvz = particles[i].vz - sun.vz;
		// Benitez and Gallardo 2008
		const double alpha = G*sun.m/(_r*_r*_r*C*C);
		const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
		const double beta = 4.*G*sun.m/_r - v2;
		const double gamma = 4.*(dvx*dx + dvy*dy + dvz*dz);

		const double dax = alpha*(beta*dx + gamma*dvx);
		const double day = alpha*(beta*dy + gamma*dvy);
		const double daz = alpha*(beta*dz + gamma*dvz);
		const double massratio = particles[i].m/particles[0].m;

		particles[i].ax += dax;
		particles[i].ay += day;
		particles[i].az += daz;
		particles[0].ax -= massratio*dax;
		particles[0].ay -= massratio*day;
		particles[0].az -= massratio*daz;
	}
}

void rebx_gr_potential(struct reb_simulation* const sim){
	// Nobili & Roxburgh 1986
	struct rebx_params_gr* rebxparams = &((struct rebx_extras*)(sim->extras))->gr;
	const double C = rebxparams->c;
	const int _N_real = sim->N - sim->N_var;
	const double G = sim->G;
	struct reb_particle* const particles = sim->particles;
	
	const struct reb_particle sun = particles[0];
	const double prefac1 = 6.*(G*sun.m)*(G*sun.m)/(C*C);
	for (int i=1; i<_N_real; i++){
		const double dx = particles[i].x - sun.x;
		const double dy = particles[i].y - sun.y;
		const double dz = particles[i].z - sun.z;
		const double r2 = dx*dx + dy*dy + dz*dz;
		const double prefac = prefac1/(r2*r2);
		
		particles[i].ax -= prefac*dx;
		particles[i].ay -= prefac*dy;
		particles[i].az -= prefac*dz;
	}
}

