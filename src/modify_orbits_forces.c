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
#include <stdlib.h>
#include "reboundx.h"

void rebx_modify_orbits_forces(struct reb_simulation* const sim){
	struct rebx_extras* rebx = sim->extras;
	struct rebx_params_modify_orbits modparams = rebx->modify_orbits_forces;

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
		fprintf(stderr, "Coordinate system not set in modify_orbits_forces?! \n");
		exit(1);
	}

	for(int i=sim->N-1;i>0;--i){
		struct reb_particle* p = &(sim->particles[i]);
		const double tau_a = rebx_get_tau_a(p);
		const double tau_e = rebx_get_tau_e(p);
		const double tau_inc = rebx_get_tau_inc(p);
		//const double tau_omega = rebx_get_tau_omega(p);
		//const double tau_Omega = rebx_get_tau_Omega(p);

		const double dvx = p->vx - com.vx;
		const double dvy = p->vy - com.vy;
		const double dvz = p->vz - com.vz;

		p->ax +=  dvx/(2.*tau_a);
		p->ay +=  dvy/(2.*tau_a);
		p->az +=  dvz/(2.*tau_a);


		if (tau_e < INFINITY || tau_inc < INFINITY){// || diskmass != 0.){ 	// need h and e vectors for both types
			const double dx = p->x-com.x;
			const double dy = p->y-com.y;
			const double dz = p->z-com.z;
			const double r = sqrt ( dx*dx + dy*dy + dz*dz );
			const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
			const double prefac = 2*vr/r;
			p->ax += prefac*dx/tau_e;
			p->ay += prefac*dy/tau_e;
			p->az += prefac*dz/tau_e + 2.*dvz/tau_inc;

			/*Alternate damping schemes:
			 * const double mu = sim->G*(com.m + p->m);
			 * const double hx = dy*dvz - dz*dvy;
			 * const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
			const double hy = dz*dvx - dx*dvz;
			const double hz = dx*dvy - dy*dvx;
			const double h = sqrt ( hx*hx + hy*hy + hz*hz );
			const double ex = 1./mu*( (v*v-mu/r)*dx - r*vr*dvx );
			const double ey = 1./mu*( (v*v-mu/r)*dy - r*vr*dvy );
			const double ez = 1./mu*( (v*v-mu/r)*dz - r*vr*dvz );
			const double e = sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity
				 *
				 * const double a = -mu/( v*v - 2.*mu/r ); // semi major axis
				const double prefac1 = 1./tau_e/1.5*(1.+modparams->e_damping_p/2.*e*e);
				const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))/tau_e/1.5;

				p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
				p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
				p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;*/



				/*const double prefac = -(hx*hx + hy*hy)/h/h/tau_inc;
				p->ax += prefac*dvx;
				p->ay += prefac*dvy;
				p->az += prefac*dvz;*/

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
		if(modparams.coordinates == JACOBI){
			rebxtools_update_com_without_particle(&com, p);
		}
	}
}
