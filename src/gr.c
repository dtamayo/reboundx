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
    const unsigned int _gravity_ignore_10 = sim->gravity_ignore_10;
    
    double aoverm10x, aoverm10y, aoverm10z;

    if (_gravity_ignore_10){
        const double dx = particles[0].x - particles[1].x;
        const double dy = particles[0].y - particles[1].y;
        const double dz = particles[0].z - particles[1].z;
        const double softening2 = sim->softening*sim->softening;
        const double r2 = dx*dx + dy*dy + dz*dz + softening2;
        const double r = sqrt(r2);
        const double prefac = G/(r2*r);
        
        aoverm10x = prefac*dx;
        aoverm10y = prefac*dy;
        aoverm10z = prefac*dz;
    }
	
    const struct reb_particle source = particles[source_index];
    for (int i=0; i<_N_real; i++){
        if(i == source_index){
            continue;
        }
		const struct reb_particle pi = particles[i];
		
		const double dx = pi.x - source.x;
		const double dy = pi.y - source.y;
		const double dz = pi.z - source.z;
		const double r2 = dx*dx + dy*dy + dz*dz;
		const double r = sqrt(r2);
		const double vx = pi.vx;
		const double vy = pi.vy;
		const double vz = pi.vz;
		const double v2 = vx*vx + vy*vy + vz*vz;
        double ax = pi.ax;
        double ay = pi.ay;
        double az = pi.az;
        if(_gravity_ignore_10 && i==1){
            ax += particles[0].m*aoverm10x;
            ay += particles[0].m*aoverm10y;
            az += particles[0].m*aoverm10z;
        }
        if(_gravity_ignore_10 && i==0){
            ax -= particles[1].m*aoverm10x;
            ay -= particles[1].m*aoverm10y;
            az -= particles[1].m*aoverm10z;
        }
        
		const double a1_x = (mu*mu*dx/(r2*r2) - 3.*mu*v2*dx/(2.*r2*r))/(C*C);
		const double a1_y = (mu*mu*dy/(r2*r2) - 3.*mu*v2*dy/(2.*r2*r))/(C*C);
		const double a1_z = (mu*mu*dz/(r2*r2) - 3.*mu*v2*dz/(2.*r2*r))/(C*C);

		const double va = vx*ax + vy*ay + vz*az;
		const double rv = dx*vx + dy*vy + dz*vz;
	
        ax = a1_x-(va*vx + v2*ax/2. + 3.*mu*(ax*r-vx*rv/r)/r2)/(C*C);
        ay = a1_y-(va*vy + v2*ay/2. + 3.*mu*(ay*r-vy*rv/r)/r2)/(C*C);
        az = a1_z-(va*vz + v2*az/2. + 3.*mu*(az*r-vz*rv/r)/r2)/(C*C);

		particles[i].ax += ax;
		particles[i].ay += ay; 
		particles[i].az += az;
		particles[source_index].ax -= pi.m/source.m*ax; 
		particles[source_index].ay -= pi.m/source.m*ay; 
		particles[source_index].az -= pi.m/source.m*az; 

    }	
}

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
	const struct reb_particle source = particles[source_index];
	for (int i=0;i<_N_real;i++){
		struct reb_particle pi = particles[i];
		if (i != source_index){
			double dx = pi.x - source.x;
			double dy = pi.y - source.y;
			double dz = pi.z - source.z;
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
