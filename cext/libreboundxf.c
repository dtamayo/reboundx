#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "libreboundxf.h"
#include "xftools.h"
#include "rebound.h"

//disk parameters for precession
double gam;
double Rc;
double diskmass;
double alpha_over_rGM0;
double podot; // pericenter precession at r = Rc

// pointers for damping timescales
double *tau_a = NULL;
double *tau_e = NULL;
double *tau_i = NULL;
double *tau_po = NULL;

static double e_damping_p = 0; // p parameter from Deck & Batygin (2015) for how e-damping
// is coupled to a-damping at order e^2
// p = 1 : e-damping at const angular momentum.  p = 0 : no contribution to a-damping
// equal to p/3 with p defined as in Goldreich & Schlichting 2014

double get_e_damping_p(void){
	return e_damping_p;
}

void set_e_damping_p(double val){
	e_damping_p = val;
}

void forces(struct reb_simulation* const sim){
	struct reb_particle com = sim->particles[0]; // calculate add. forces w.r.t. center of mass
	for(int i=1;i<sim->N;i++){
		struct reb_particle* p = &(sim->particles[i]);
		const double dvx = p->vx - com.vx;
		const double dvy = p->vy - com.vy;
		const double dvz = p->vz - com.vz;

		if (tau_a[i] != 0.){
			p->ax -=  dvx/(2.*tau_a[i]);
			p->ay -=  dvy/(2.*tau_a[i]);
			p->az -=  dvz/(2.*tau_a[i]);
		}

		if (tau_e[i] != 0. || tau_i[i]!= 0. || diskmass != 0.){ 	// need h and e vectors for both types
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
			if (tau_e[i] != 0.){	// Eccentricity damping
				/*const double a = -mu/( v*v - 2.*mu/r );			// semi major axis
				const double prefac1 = 1./tau_e[i]/1.5*(1.+e_damping_p/2.*e*e);
				const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))/tau_e[i]/1.5;*/

				p->ax += -2/tau_e[i]*vr*dx/r;
				p->ay += -2/tau_e[i]*vr*dy/r;
				p->az += -2/tau_e[i]*vr*dz/r;
				/*p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
				p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
				p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;*/
			}
			if (tau_i[i]!=0){		// Inclination damping
				p->az += -2.*dvz/tau_i[i];
				const double prefac = (hx*hx + hy*hy)/h/h/tau_i[i];
				p->ax += prefac*dvx;
				p->ay += prefac*dvy;
				p->az += prefac*dvz;
			}
			if (diskmass != 0.) {
				double a_over_r = -sim->G*sim->particles[0].m*alpha_over_rGM0*pow(Rc/r,gam)/r + sim->G*diskmass/r/r/r; 	// radial disk force after removing piece from adding the disk into the sun
				p->ax += a_over_r*dx;									// rhat has components x/r xhat + y/r yhat + z/r zhat
				p->ay += a_over_r*dy;
				p->az += a_over_r*dz;

				sim->particles[0].ax -= p->m/sim->particles[0].m*a_over_r*dx;		// add back reactions onto the star (if forces are equal, accelerations differ by -mass ratio)
				sim->particles[0].ay -= p->m/sim->particles[0].m*a_over_r*dy;
				sim->particles[0].az -= p->m/sim->particles[0].m*a_over_r*dz;
			}
		}
		com = xftools_get_com(com,sim->particles[i]);
	}
	xftools_move_to_com(sim->particles, sim->N);
}

void modify_elements(struct reb_simulation* const sim){
	struct reb_particle com = sim->particles[0];
	for(int i=1;i<sim->N;i++){
		struct reb_particle *p = &(sim->particles[i]);
		struct reb_orbit o = xftools_p2orbit(sim->G, sim->particles[i], com);
	    double da = 0.;
		double de = 0.;
		double dpo = 0.;	
		if (tau_a[i] != 0.){
			da += -o.a*sim->dt/tau_a[i]; 
		}
		
		if (tau_e[i] != 0.){
			de += -o.e*sim->dt/tau_e[i];
			da += -2.*o.a*o.e*o.e*e_damping_p*sim->dt/tau_e[i];
		}

		if (tau_po[i] != 0.){
			dpo += 2*M_PI*sim->dt/tau_po[i];
		}

		o.a += da;
		o.e += de;
		o.omega += dpo;

		xftools_orbit2p(&sim->particles[i], sim->G, &com, o); 

		com = xftools_get_com(com, sim->particles[i]);
	}
}

static void xf_init(int N){ // only used internally
	if(tau_a == NULL){	tau_a = calloc(sizeof(double),N);}
	if(tau_e == NULL){	tau_e = calloc(sizeof(double),N);}
	if(tau_i == NULL){	tau_i = calloc(sizeof(double),N);}
	if(tau_po == NULL){  tau_po = calloc(sizeof(double),N);}
}
	
void set_migration(double *_tau_a, int N){
	/*if(N > 0 && N != N){
		printf("A previous call to reboundxf used a different number of particles, which is not supported in the current implementation.  Please improve me!\n");
		exit(1);
	}*/
	
	//N = N;

	xf_init(N);
	for(int i=0; i<N; ++i){
		tau_a[i] = _tau_a[i];
	}
}

void set_e_damping(double *_tau_e, int N){
	/*if(N > 0 && N != N){
		printf("A previous call to reboundxf used a different number of particles, which is not supported in the current implementation.  Please improve me!\n");
		exit(1);
	}*/
	
	//N = N;
	xf_init(N);
	for(int i=0; i<N; ++i){
		tau_e[i] = _tau_e[i];
	}
}

void set_i_damping(double *_tau_i, int N){
	/*if(N > 0 && N != N){
		printf("A previous call to reboundxf used a different number of particles, which is not supported in the current implementation.  Please improve me!\n");
		exit(1);
	}*/
	
	//N = N;
	xf_init(N);
	if(tau_i == NULL){	tau_i = calloc(sizeof(double),N);}
	for(int i=0; i<N; ++i){
		tau_i[i] = _tau_i[i];
	}
}

void set_peri_precession(double *_tau_po, int N){
	xf_init(N);
	for(int i=0; i<N; ++i){
		tau_po[i] = _tau_po[i];
	}
}

/* not yet implemented
void set_peri_precession(double _gam, double _Rc, double _podot, int N){
	xf_init(N);
	Rc = _Rc;
	gam = _gam;
	podot = _podot; // as a fraction of the mean motion

	alpha_over_rGM0 = Rc/Rc*podot/(1.-gam/2.);
	diskmass = 3.65557*sim->particles[0].m*podot;
	
	sim->particles[0].m += diskmass;
}
*/
void reset(){
	free(tau_a);
	tau_a = NULL;
	free(tau_e);
	tau_e = NULL;
	free(tau_i);
	tau_i = NULL;
	Rc=0.;
	gam = 0.;
	podot = 0.;
	alpha_over_rGM0 = 0.;
	diskmass = 0.;

}

