#include <stdlib.h>
#include <string.h>
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

void rebxf_check_N(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	if(xf->allocatedN < sim->N) {
		xf->allocatedN = sim->N;
		xf->tau_a = realloc(xf->tau_a, sizeof(double)*xf->allocatedN);
		xf->tau_e = realloc(xf->tau_e, sizeof(double)*xf->allocatedN);
		xf->tau_inc = realloc(xf->tau_inc, sizeof(double)*xf->allocatedN);
		xf->tau_pomega = realloc(xf->tau_pomega, sizeof(double)*xf->allocatedN);
	}
}

double* rebxf_get_tau_a(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	return xf->tau_a;	
}

void rebxf_set_tau_a(struct reb_simulation* sim, double* tau_a){
	rebxf_check_N(sim);
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	memcpy(xf->tau_a, tau_a, sizeof(double)*sim->N);	
}

double* rebxf_get_tau_e(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	return xf->tau_e;	
}

void rebxf_set_tau_e(struct reb_simulation* sim, double* tau_e){
	rebxf_check_N(sim);
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	memcpy(xf->tau_e, tau_e, sizeof(double)*sim->N);	
}

double* rebxf_get_tau_inc(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	return xf->tau_inc;	
}

void rebxf_set_tau_inc(struct reb_simulation* sim, double* tau_inc){
	rebxf_check_N(sim);
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	memcpy(xf->tau_inc, tau_inc, sizeof(double)*sim->N);	
}

double* rebxf_get_tau_pomega(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	return xf->tau_pomega;	
}

void rebxf_set_tau_pomega(struct reb_simulation* sim, double* tau_pomega){
	rebxf_check_N(sim);
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	memcpy(xf->tau_pomega, tau_pomega, sizeof(double)*sim->N);	
}

struct rebxf_params* rebxf_addxf(struct reb_simulation* sim){
	struct rebxf_params* xf = (struct rebxf_params*) malloc(sizeof(struct rebxf_params));
	xf->tau_a = NULL;
	xf->tau_e = NULL;
	xf->tau_inc = NULL;
	xf->tau_pomega = NULL;

	xf->e_damping_p = 0.;

	xf->allocatedN = sim->N;
	xf->tau_a = calloc(xf->allocatedN, sizeof(double));
	xf->tau_e = calloc(xf->allocatedN, sizeof(double));
	xf->tau_inc = calloc(xf->allocatedN, sizeof(double));
	xf->tau_pomega = calloc(xf->allocatedN, sizeof(double));

	sim->xf_params = (struct rebxf_params*) xf;
	return (struct rebxf_params*)sim->xf_params;
}

void rebxf_forces(struct reb_simulation* const sim){
	rebxf_check_N(sim);
	struct rebxf_params* xf = (struct rebxf_params*)xf;
	struct reb_particle com = sim->particles[0]; // calculate add. forces w.r.t. center of mass
	for(int i=1;i<sim->N;i++){
		struct reb_particle* p = &(sim->particles[i]);
		const double dvx = p->vx - com.vx;
		const double dvy = p->vy - com.vy;
		const double dvz = p->vz - com.vz;

		if (xf->tau_a[i] != 0.){
			p->ax -=  dvx/(2.*xf->tau_a[i]);
			p->ay -=  dvy/(2.*xf->tau_a[i]);
			p->az -=  dvz/(2.*xf->tau_a[i]);
		}

		if (xf->tau_e[i] != 0. || xf->tau_inc[i]!= 0.){// || diskmass != 0.){ 	// need h and e vectors for both types
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
			if (xf->tau_inc[i]!=0){		// Inclination damping
				p->az += -2.*dvz/xf->tau_inc[i];
				const double prefac = (hx*hx + hy*hy)/h/h/xf->tau_inc[i];
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
		com = xftools_get_com(com,sim->particles[i]);
	}
	xftools_move_to_com(sim->particles, sim->N);
}

void rebxf_modify_elements(struct reb_simulation* const sim){
	rebxf_check_N(sim);
	struct rebxf_params* xf = (struct rebxf_params*)sim->xf_params;
	struct reb_particle com = sim->particles[0];
	for(int i=1;i<sim->N;i++){
		struct reb_particle *p = &(sim->particles[i]);
		struct reb_orbit o = xftools_p2orbit(sim->G, sim->particles[i], com);
	    double da = 0.;
		double de = 0.;
		double dpo = 0.;	
		if (xf->tau_a[i] != 0.){
			da += -o.a*sim->dt/xf->tau_a[i]; 
		}
		
		if (xf->tau_e[i] != 0.){
			de += -o.e*sim->dt/xf->tau_e[i];
			da += -2.*o.a*o.e*o.e*xf->e_damping_p*sim->dt/xf->tau_e[i];
		}

		if (xf->tau_pomega[i] != 0.){
			dpo += 2*M_PI*sim->dt/xf->tau_pomega[i]*(1.+sin(o.omega));
		}

		o.a += da;
		o.e += de;
		o.omega += dpo;

		xftools_orbit2p(&sim->particles[i], sim->G, &com, o); 
		com = xftools_get_com(com, sim->particles[i]);
	}
	xftools_move_to_com(sim->particles, sim->N);
}
