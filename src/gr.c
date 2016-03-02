/**
 * @file    gr.c
 * @brief   Post-newtonian general relativity corrections arising from a single massive body
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
#include "gr.h"
#include "rebound.h"
#include "reboundx.h"

struct rebx_params_gr* rebx_add_gr(struct rebx_extras* rebx, int source_index, double c){
	struct rebx_params_gr* params = malloc(sizeof(*params));
	params->c = c;
    params->source_index = source_index;
    int force_is_velocity_dependent = 1;
    rebx_add_force(rebx, params, "gr", rebx_gr, force_is_velocity_dependent);
    return params;
}

void rebx_gr(struct reb_simulation* const sim, struct rebx_effect* gr){
    const struct rebx_params_gr* const params = gr->paramsPtr;
    const double C = params->c;
    const int source_index = params->source_index;
    const int _N_real = sim->N - sim->N_var;
    const double G = sim->G;
    struct reb_particle* const particles = sim->particles;
	const double mu = G*particles[source_index].m;
    for (int i=0; i<_N_real; i++){
        if(i == source_index){
            continue;
        }
		struct reb_particle pi = particles[i];
		
		double dx = pi.x;
		double dy = pi.y;
		double dz = pi.z;
		double r2 = dx*dx + dy*dy + dz*dz;
		double r = sqrt(r2);
		double vx = pi.vx;
		double vy = pi.vy;
		double vz = pi.vz;
		double v2 = vx*vx + vy*vy + vz*vz;

		double a1_x = (mu*mu*dx/(r2*r2) - 3.*mu*v2*dx/(2.*r2*r))/(C*C);
		double a1_y = (mu*mu*dy/(r2*r2) - 3.*mu*v2*dy/(2.*r2*r))/(C*C);
		double a1_z = (mu*mu*dz/(r2*r2) - 3.*mu*v2*dz/(2.*r2*r))/(C*C);

		double va = vx*pi.ax + vy*pi.ay + vz*pi.az;
		double rv = dx*vx + dy*vy + dz*vz;
		
		particles[i].ax += a1_x-(va*vx + v2*pi.ax/2. + 3.*mu*(pi.ax*r-vx*rv/r)/r2)/(C*C);
		particles[i].ay += a1_y-(va*vy + v2*pi.ay/2. + 3.*mu*(pi.ay*r-vy*rv/r)/r2)/(C*C);
		particles[i].az += a1_z-(va*vz + v2*pi.az/2. + 3.*mu*(pi.az*r-vz*rv/r)/r2)/(C*C);
    }	
}
/*void rebx_gr(struct reb_simulation* const sim, struct rebx_effect* gr){
    const struct rebx_params_gr* const params = gr->paramsPtr;
    const double C = params->c;
    const int source_index = params->source_index;
    const int _N_real = sim->N - sim->N_var;
    const double G = sim->G;
    struct reb_particle* const particles = sim->particles;
    const struct reb_particle source = sim->particles[source_index];
    for (int i=0; i<_N_real; i++){
        if(i == source_index){
            continue;
        }
        const struct reb_particle p = sim->particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double _r = sqrt(r2);
        const double dvx = p.vx - source.vx;
        const double dvy = p.vy - source.vy;
        const double dvz = p.vz - source.vz;
        // Benitez and Gallardo 2008
        const double alpha = G*source.m/(_r*_r*_r*C*C);
        const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
        const double beta = 4.*G*source.m/_r - v2;
        const double gamma = 4.*(dvx*dx + dvy*dy + dvz*dz);

        const double dax = alpha*(beta*dx + gamma*dvx);
        const double day = alpha*(beta*dy + gamma*dvy);
        const double daz = alpha*(beta*dz + gamma*dvz);
        //const double massratio = particles[i].m/source.m;

        particles[i].ax += dax;
        particles[i].ay += day;
        particles[i].az += daz;
        //particles[source_index].ax -= massratio*dax;
        //particles[source_index].ay -= massratio*day;
        //particles[source_index].az -= massratio*daz;
    }
}
*/

double rebx_gr_hamiltonian(const struct reb_simulation* const sim, const struct rebx_params_gr* const params){ 
    const double C = params->c;
    const int source_index = params->source_index;
	const struct reb_particle* const particles = sim->particles;
	const int _N_real = sim->N - sim->N_var;
	const double G = sim->G;
	const double mu = G*particles[source_index].m;

	double e_kin = 0.;
	double e_pot = 0.;
	double e_pn  = 0.;
	struct reb_particle source = particles[source_index];
	for (int i=0;i<_N_real;i++){
		struct reb_particle pi = particles[i];
		if (i != source_index){
			double dx = pi.x;
			double dy = pi.y;
			double dz = pi.z;
			double r2 = dx*dx + dy*dy + dz*dz;
			double r = sqrt(r2);

			double vx = pi.vx;
			double vy = pi.vy;
			double vz = pi.vz;
			double v2 = vx*vx + vy*vy + vz*vz;
			
			double A = 1. - (v2/2. + 3.*mu/r)/(C*C);
			double B = sqrt(v2)/A;
			double v_tilde2 = B*B;

			e_kin += 0.5*pi.m*v_tilde2;
			e_pn += (mu*mu*pi.m/(2.*r2) - v_tilde2*v_tilde2*pi.m/8. - 3.*mu*v_tilde2*pi.m/(2.*r))/(C*C);
		}		
		else{
			double source_v2 = source.vx*source.vx + source.vy*source.vy + source.vz*source.vz;
			e_kin += 0.5 * source.m * source_v2;
		}
		for (int j=i+1; j<_N_real; j++){
			struct reb_particle pj = particles[j];
			double dx = pi.x - pj.x;
			double dy = pi.y - pj.y;
			double dz = pi.z - pj.z;	
			double r = sqrt(dx*dx + dy*dy + dz*dz);

			e_pot -= G*pi.m*pj.m/r;
		}
	}
	return e_kin + e_pot + e_pn;
}
