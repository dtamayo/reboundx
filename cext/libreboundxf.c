#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "libreboundxf.h"
#include "xftools.h"

//disk parameters for precession
/*double gam;
double Rc;
double diskmass;
double alpha_over_rGM0;
double podot; // pericenter precession at r = Rc
*/
// pointers for damping timescales

struct rebxf_params* rebxf_addxf(struct reb_simulation* const sim){
	struct rebxf_params* xf = calloc(1, sizeof(struct rebxf_params));
	xf->tau_a = NULL;
	xf->tau_e = NULL;
	xf->tau_i = NULL;
	xf->tau_po = NULL;

	xf->e_damping_p = 0.;
	xf->sim = sim;

	xf->allocatedN = 128;
	xf->tau_a = calloc(xf->allocatedN, sizeof(double));
	xf->tau_e = calloc(xf->allocatedN, sizeof(double));
	xf->tau_i = calloc(xf->allocatedN, sizeof(double));
	xf->tau_po = calloc(xf->allocatedN, sizeof(double));

	return xf;
}

inline void rebxf_checkN(struct rebxf_params* const xf){
	if(xf->allocatedN <= xf->sim->N) {
		xf->allocatedN += 128;
		xf->tau_a = realloc(xf->tau_a, sizeof(double)*xf->allocatedN);
		xf->tau_e = realloc(xf->tau_e, sizeof(double)*xf->allocatedN);
		xf->tau_i = realloc(xf->tau_i, sizeof(double)*xf->allocatedN);
		xf->tau_po = realloc(xf->tau_po, sizeof(double)*xf->allocatedN);
	}
}

void rebxf_forces(struct rebxf_params* const xf){
	struct reb_particle com = xf->sim->particles[0]; // calculate add. forces w.r.t. center of mass
	for(int i=1;i<xf->sim->N;i++){
		struct reb_particle* p = &(xf->sim->particles[i]);
		const double dvx = p->vx - com.vx;
		const double dvy = p->vy - com.vy;
		const double dvz = p->vz - com.vz;

		if (xf->tau_a[i] != 0.){
			p->ax -=  dvx/(2.*xf->tau_a[i]);
			p->ay -=  dvy/(2.*xf->tau_a[i]);
			p->az -=  dvz/(2.*xf->tau_a[i]);
		}

		if (xf->tau_e[i] != 0. || xf->tau_i[i]!= 0.){// || diskmass != 0.){ 	// need h and e vectors for both types
			const double mu = xf->sim->G*(com.m + p->m);
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
			//printf("%.14f\t%.2e\n", vr/v, e);
			if (xf->tau_e[i] != 0.){	// Eccentricity damping
				/*const double a = -mu/( v*v - 2.*mu/r );			// semi major axis
				const double prefac1 = 1./xf->tau_e[i]/1.5*(1.+e_damping_p/2.*e*e);
				const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))/xf->tau_e[i]/1.5;*/

				p->ax += -2/xf->tau_e[i]*vr*dx/r;
				p->ay += -2/xf->tau_e[i]*vr*dy/r;
				p->az += -2/xf->tau_e[i]*vr*dz/r;
				/*p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
				p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
				p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;*/
			}
			if (xf->tau_i[i]!=0){		// Inclination damping
				p->az += -2.*dvz/xf->tau_i[i];
				const double prefac = (hx*hx + hy*hy)/h/h/xf->tau_i[i];
				p->ax += prefac*dvx;
				p->ay += prefac*dvy;
				p->az += prefac*dvz;
			}
			/*if (diskmass != 0.) {
				double a_over_r = -xf->sim->G*xf->sim->particles[0].m*alpha_over_rGM0*pow(Rc/r,gam)/r + xf->sim->G*diskmass/r/r/r; 	// radial disk force after removing piece from adding the disk into the sun
				p->ax += a_over_r*dx;									// rhat has components x/r xhat + y/r yhat + z/r zhat
				p->ay += a_over_r*dy;
				p->az += a_over_r*dz;

				xf->sim->particles[0].ax -= p->m/xf->sim->particles[0].m*a_over_r*dx;		// add back reactions onto the star (if forces are equal, accelerations differ by -mass ratio)
				xf->sim->particles[0].ay -= p->m/xf->sim->particles[0].m*a_over_r*dy;
				xf->sim->particles[0].az -= p->m/xf->sim->particles[0].m*a_over_r*dz;
			}*/
		}
		com = xftools_get_com(com,xf->sim->particles[i]);
	}
	xftools_move_to_com(xf->sim->particles, xf->sim->N);
}

void rebxf_modify_elements(struct rebxf_params* const xf){
	struct reb_particle com = xf->sim->particles[0];
	for(int i=1;i<xf->sim->N;i++){
		struct reb_particle *p = &(xf->sim->particles[i]);
		struct reb_orbit o = xftools_p2orbit(xf->sim->G, xf->sim->particles[i], com);
	    double da = 0.;
		double de = 0.;
		double dpo = 0.;	
		if (xf->tau_a[i] != 0.){
			da += -o.a*xf->sim->dt/xf->tau_a[i]; 
		}
		
		if (xf->tau_e[i] != 0.){
			de += -o.e*xf->sim->dt/xf->tau_e[i];
			da += -2.*o.a*o.e*o.e*xf->e_damping_p*xf->sim->dt/xf->tau_e[i];
		}

		if (xf->tau_po[i] != 0.){
			dpo += 2*M_PI*xf->sim->dt/xf->tau_po[i]*(1.+sin(o.omega));
		}

		o.a += da;
		o.e += de;
		o.omega += dpo;

		xftools_orbit2p(&xf->sim->particles[i], xf->sim->G, &com, o); 
		com = xftools_get_com(com, xf->sim->particles[i]);
	}
	xftools_move_to_com(xf->sim->particles, xf->sim->N);
}
