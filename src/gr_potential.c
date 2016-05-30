/**
 * @file    gr_potential.c
 * @brief   Post-newtonian general relativity corrections using a simple potential that gets the pericenter precession right.
 * @author  Pengshuai (Sam) Shi, Hanno Rein, Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Pengshuai (Sam) Shi, Hanno Rein, Dan Tamayo
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gr_potential.h"
#include "rebound.h"
#include "reboundx.h"

struct rebx_params_gr_potential* rebx_add_gr_potential(struct rebx_extras* rebx, int source_index, double c){
	struct rebx_params_gr_potential* params = malloc(sizeof(*params));
	params->c = c;
    params->source_index = source_index;
    int force_is_velocity_dependent = 0;
    rebx_add_force(rebx, params, "gr_potential", rebx_gr_potential, force_is_velocity_dependent);
    return params;
}

void rebx_gr_potential(struct reb_simulation* const sim, struct rebx_effect* gr){
    // Nobili & Roxburgh 1986
    const struct rebx_params_gr_potential* const params = gr->paramsPtr;
    const double C = params->c;
    const int source_index = params->source_index;
    const int _N_real = sim->N - sim->N_var;
    const double G = sim->G;
    struct reb_particle* const particles = sim->particles;
    const struct reb_particle source = sim->particles[source_index];
    
    const double prefac1 = 6.*(G*source.m)*(G*source.m)/(C*C);
    for (int i=0; i<_N_real; i++){
        if(i == source_index){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double prefac = prefac1/(r2*r2);
        
        particles[i].ax -= prefac*dx;
        particles[i].ay -= prefac*dy;
        particles[i].az -= prefac*dz;
        particles[source_index].ax += p.m/source.m*prefac*dx;
        particles[source_index].ay += p.m/source.m*prefac*dy;
        particles[source_index].az += p.m/source.m*prefac*dz;
    }
}

double rebx_gr_potential_hamiltonian(const struct reb_simulation* const sim, const struct rebx_params_gr_potential* const params){ 
    const double C = params->c;
    const int source_index = params->source_index;
	const struct reb_particle* const particles = sim->particles;
	const int _N_real = sim->N - sim->N_var;
	const double G = sim->G;
    const struct reb_particle source = particles[source_index];
	const double mu = G*source.m;
    const double prefac = 3.*mu*mu/(C*C);
    double H = reb_tools_energy(sim);

	for (int i=0;i<_N_real;i++){
		if(i == source_index){
            continue;
        }
        struct reb_particle pi = particles[i];
        double dx = pi.x - source.x;
        double dy = pi.y - source.y;
        double dz = pi.z - source.z;
        double r2 = dx*dx + dy*dy + dz*dz;
        H -= prefac*pi.m/r2;
    }		
	
    return H;
}

