#include <stdlib.h>
#include <math.h>
#include <stdio.h>

struct particle {
	double x;	/**< x-position of the particle. */
	double y;	/**< y-position of the particle. */
	double z;	/**< z-position of the particle. */
	double vx;	/**< x-velocity of the particle. */
	double vy;	/**< y-velocity of the particle. */
	double vz;	/**< z-velocity of the particle. */
	double ax;	/**< x-acceleration of the particle. */
	double ay;	/**< y-acceleration of the particle. */
	double az;	/**< z-acceleration of the particle. */
	double m;	/**< Mass of the particle. */
#ifndef COLLISIONS_NONE
	double r; 	/**< Radius of the particle. */
	double lastcollision;	/**< Last time the particle had a physical collision. */
#endif // COLLISIONS_NONE
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
	struct cell* c;		/**< Pointer to the cell the particle is currently in. */
#endif // TREE
};

//static int N = 0; // private version to make sure we don't overwrite/change what's in simulation

//disk parameters for precession
double gam;
double Rc;
double diskmass;
double alpha_over_GM0;
double podot; // pericenter precession at r = Rc

// pointers for damping timescales
double *tau_a = NULL;
double *tau_e = NULL;
double *tau_i = NULL;

//not yet implemented
//double e_damping_p; // p parameter from Goldreich & Schlichting 2014 for how e-damping
// contributes to a-damping at order e^2
// p = 3 : e-damping at const angular momentum.  p = 0 : no contribution to a-damping


void move_to_com(struct particle* particles, int _N){
	double m = 0;
	double x = 0;
	double y = 0;
	double z = 0;
	double vx = 0;
	double vy = 0;
	double vz = 0;
	for (int i=0;i<_N;i++){
		struct particle p = particles[i];
		m  += p.m;
		x  += p.x*p.m;
		y  += p.y*p.m;
		z  += p.z*p.m;
		vx += p.vx*p.m;
		vy += p.vy*p.m;
		vz += p.vz*p.m;
	}
	x /= m;
	y /= m;
	z /= m;
	vx /= m;
	vy /= m;
	vz /= m;
	for (int i=0;i<_N;i++){
		particles[i].x  -= x;
		particles[i].y  -= y;
		particles[i].z  -= z;
		particles[i].vx -= vx;
		particles[i].vy -= vy;
		particles[i].vz -= vz;
	}
}

struct particle get_com(struct particle p1, struct particle p2){
	p1.x   = p1.x*p1.m + p2.x*p2.m;		
	p1.y   = p1.y*p1.m + p2.y*p2.m;
	p1.z   = p1.z*p1.m + p2.z*p2.m;
	p1.vx  = p1.vx*p1.m + p2.vx*p2.m;
	p1.vy  = p1.vy*p1.m + p2.vy*p2.m;
	p1.vz  = p1.vz*p1.m + p2.vz*p2.m;
	p1.m  += p2.m;
	if (p1.m>0.){
		p1.x  /= p1.m;
		p1.y  /= p1.m;
		p1.z  /= p1.m;
		p1.vx /= p1.m;
		p1.vy /= p1.m;
		p1.vz /= p1.m;
	}
	return p1;
}

void disk_forces(struct particle* particles, double t, double dt, double G, int _N, int N_megno){	
	printf("%d\n", _N);
	struct particle com = particles[0]; // calculate add. forces w.r.t. center of mass
	for(int i=1;i<_N;i++){
		struct particle* p = &(particles[i]);
		const double dvx = p->vx - com.vx;
		const double dvy = p->vy - com.vy;
		const double dvz = p->vz - com.vz;

		if (tau_a[i] != 0.){
			p->ax -=  dvx/(2.*tau_a[i]);
			p->ay -=  dvy/(2.*tau_a[i]);
			p->az -=  dvz/(2.*tau_a[i]);
		}

		if (tau_e[i] != 0. || tau_i[i]!= 0. || diskmass != 0.){ 	// need h and e vectors for both types
			const double mu = G*(com.m + p->m);
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
			//const double e = sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity

			if (tau_e[i] != 0.){	// Eccentricity damping
				//const double a = -mu/( v*v - 2.*mu/r );			// semi major axis
				//const double prefac1 = 1./tau_e[i]/1.5*(1.+e_damping_p/2.*e*e);
				//const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))/tau_e[i]/1.5;

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
				double a_over_r = -G*particles[0].m*alpha_over_GM0*pow(Rc/r,gam)/r + G*diskmass/r/r/r; 	// radial disk force after removing piece from adding the disk into the sun
				p->ax += a_over_r*dx;									// rhat has components x/r xhat + y/r yhat + z/r zhat
				p->ay += a_over_r*dy;
				p->az += a_over_r*dz;

				particles[0].ax -= p->m/particles[0].m*a_over_r*dx;		// add back reactions onto the star (if forces are equal, accelerations differ by -mass ratio)
				particles[0].ay -= p->m/particles[0].m*a_over_r*dy;
				particles[0].az -= p->m/particles[0].m*a_over_r*dz;
			}
		}
		com = get_com(com,particles[i]);
	}
	move_to_com(particles, _N);
}

void set_migration(double *_tau_a, int _N){
	/*if(N > 0 && N != _N){
		printf("A previous call to reboundxf used a different number of particles, which is not supported in the current implementation.  Please improve me!\n");
		exit(1);
	}*/
	
	//N = _N;
	if(tau_a == NULL){	tau_a = calloc(sizeof(double),_N);}
	
	for(int i=0; i<_N; ++i){
		tau_a[i] = _tau_a[i];
	}
}

void set_e_damping(double *_tau_e, int _N){
	/*if(N > 0 && N != _N){
		printf("A previous call to reboundxf used a different number of particles, which is not supported in the current implementation.  Please improve me!\n");
		exit(1);
	}*/
	
	//N = _N;
	if(tau_e == NULL){	tau_e = calloc(sizeof(double),_N);}
	for(int i=0; i<_N; ++i){
		tau_e[i] = _tau_e[i];
	}
}

void set_i_damping(double *_tau_i, int _N){
	/*if(N > 0 && N != _N){
		printf("A previous call to reboundxf used a different number of particles, which is not supported in the current implementation.  Please improve me!\n");
		exit(1);
	}*/
	
	//N = _N;
	if(tau_i == NULL){	tau_i = calloc(sizeof(double),_N);}
	for(int i=0; i<_N; ++i){
		tau_i[i] = _tau_i[i];
	}
}
/* not yet implemented
void set_peri_precession(double _gam, double _Rc, double _podot){
	Rc = _Rc;
	gam = _gam;
	podot = _podot; // as a fraction of the mean motion

	alpha_over_GM0 = Rc/Rc*podot/(1.-gam/2.);
	diskmass = 3.65557*particles[0].m*podot;
	
	particles[0].m += diskmass;
}
*/
void reboundxf_reset(){
	free(tau_a);
	tau_a = NULL;
	free(tau_e);
	tau_e = NULL;
	free(tau_i);
	tau_i = NULL;
	Rc=0.;
	gam = 0.;
	podot = 0.;
	alpha_over_GM0 = 0.;
	diskmass = 0.;

}

