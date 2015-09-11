/**
 * @file 	modify_orbits_forces.c
 * @brief 	Add forces that orbit-average to give semimajor axis and eccentricity damping
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
#include "modify_orbits_forces.h"
#include "reboundx.h"

void rebx_modify_orbits_forces(struct reb_simulation* const sim){
	//rebx_check_N(sim);
	struct rebx_params_modify_orbits_forces* rebxparams = &((struct rebx_extras*)(sim->extras))->modify_orbits_forces;
	struct reb_particle com = sim->particles[0]; // calculate add. forces w.r.t. center of mass
	for(int i=1;i<sim->N;i++){
		struct reb_particle* p = &(sim->particles[i]);
		const double dvx = p->vx - com.vx;
		const double dvy = p->vy - com.vy;
		const double dvz = p->vz - com.vz;

		if (rebxparams->tau_a[i] != 0.){
			p->ax +=  dvx/(2.*rebxparams->tau_a[i]);
			p->ay +=  dvy/(2.*rebxparams->tau_a[i]);
			p->az +=  dvz/(2.*rebxparams->tau_a[i]);
		}

		if (rebxparams->tau_e[i] != 0. || rebxparams->tau_inc[i]!= 0.){// || diskmass != 0.){ 	// need h and e vectors for both types
			const double mu = sim->G*(com.m + p->m);
			const double dx = p->x-com.x;
			const double dy = p->y-com.y;
			const double dz = p->z-com.z;
			const double hx = dy*dvz - dz*dvy;
			const double hy = dz*dvx - dx*dvz;
			const double hz = dx*dvy - dy*dvx;
			const double h = sqrt ( hx*hx + hy*hy + hz*hz );
			const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
			const double r = sqrt ( dx*dx + dy*dy + dz*dz );
			const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
			const double ex = 1./mu*( (v*v-mu/r)*dx - r*vr*dvx );
			const double ey = 1./mu*( (v*v-mu/r)*dy - r*vr*dvy );
			const double ez = 1./mu*( (v*v-mu/r)*dz - r*vr*dvz );
			const double e = sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity
			//printf("%.14f\t%.2e\n", vr/v, rebxparams->tau_e[i]);
			if (rebxparams->tau_e[i] != 0.){	// Eccentricity damping
				/*const double a = -mu/( v*v - 2.*mu/r );			// semi major axis
				const double prefac1 = 1./rebxparams->tau_e[i]/1.5*(1.+e_damping_p/2.*e*e);
				const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))/rebxparams->tau_e[i]/1.5;*/

				p->ax += 2/rebxparams->tau_e[i]*vr*dx/r;
				p->ay += 2/rebxparams->tau_e[i]*vr*dy/r;
				p->az += 2/rebxparams->tau_e[i]*vr*dz/r;
				/*p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
				p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
				p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;*/
			}
			if (rebxparams->tau_inc[i]!=0){		// Inclination damping
				p->az -= -2.*dvz/rebxparams->tau_inc[i];
				const double prefac = -(hx*hx + hy*hy)/h/h/rebxparams->tau_inc[i];
				p->ax += prefac*dvx;
				p->ay += prefac*dvy;
				p->az += prefac*dvz;
			}
			/*if (diskmass != 0.) {
				double a_over_r = -sim->G*sim->particles[0].m*alpha_over_rGM0*pow(Rc/r,gam)/r + sim->G*diskmass/r/r/r; 	// radial disk force after removing piece from adding the disk into the sun
				p->ax += a_over_r*dx;									// rhat has components x/r xhat + y/r yhat + z/r zhat
				p->ay += a_over_r*dy;
				p->az += a_over_r*dz;

				sim->particles[0].ax -= p->m/sim->particles[0].m*a_over_r*dx;		// add back reactions onto the star (if forces are equal, accelerations differ by -mass ratio)
				sim->particles[0].ay -= p->m/sim->particles[0].m*a_over_r*dy;
				sim->particles[0].az -= p->m/sim->particles[0].m*a_over_r*dz;
			}*/
		}
		com = rebxtools_get_com_of_pair(com,sim->particles[i]);
	}
	rebxtools_move_to_com(sim);
}
