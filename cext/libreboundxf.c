#include <stdlib.h>
#include <math.h>
#include <stdio.h>

struct orbit {
	double a;
	double r;	// Radial distance from central object
	double h;	// Angular momentum
	double P;	// Orbital period
	double l;
	double e;
	double inc;
	double Omega; 	// longitude of ascending node
	double omega; 	// argument of perihelion
	double f; 	// true anomaly
};

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
/*#ifndef COLLISIONSNONE
	double r; 	
	double lastcollision;
#endif // COLLISIONSNONE*/
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
double *tau_po = NULL;

double e_damping_p; // p parameter from Goldreich & Schlichting 2014 for how e-damping
// contributes to a-damping at order e^2
// p = 3 : e-damping at const angular momentum.  p = 0 : no contribution to a-damping

struct orbit orbit_nan(){
	struct orbit o;
	o.a = NAN;
	o.r = NAN;
	o.h = NAN;
	o.P = NAN;
	o.l = NAN;
	o.e = NAN;
	o.inc = NAN;
	o.Omega = NAN;
	o.omega = NAN;
	o.f = NAN;
	return o;
}
#define MIN_REL_ERROR 1.0e-12
#define TINTY 1.E-308

struct orbit tools_p2orbit(double G, struct particle p, struct particle primary){
	struct orbit o;
	double h0,h1,h2,e0,e1,e2,n0,n1,n,er,vr,mu,ea;
	mu = G*(p.m+primary.m);
	p.x -= primary.x;
	p.y -= primary.y;
	p.z -= primary.z;
	p.vx -= primary.vx;
	p.vy -= primary.vy;
	p.vz -= primary.vz;
	h0 = (p.y*p.vz - p.z*p.vy); 			//angular momentum vector
	h1 = (p.z*p.vx - p.x*p.vz);
	h2 = (p.x*p.vy - p.y*p.vx);
	o.h = sqrt ( h0*h0 + h1*h1 + h2*h2 );		// abs value of angular moment 
	double v = sqrt ( p.vx*p.vx + p.vy*p.vy + p.vz*p.vz );
	o.r = sqrt ( p.x*p.x + p.y*p.y + p.z*p.z );
	if (o.h/(o.r*v) <= MIN_REL_ERROR){
		return orbit_nan();
	}
	vr = (p.x*p.vx + p.y*p.vy + p.z*p.vz)/o.r;
	e0 = 1./mu*( (v*v-mu/o.r)*p.x - o.r*vr*p.vx );
	e1 = 1./mu*( (v*v-mu/o.r)*p.y - o.r*vr*p.vy );
	e2 = 1./mu*( (v*v-mu/o.r)*p.z - o.r*vr*p.vz );
 	o.e = sqrt( e0*e0 + e1*e1 + e2*e2 );		// eccentricity
	o.a = -mu/( v*v - 2.*mu/o.r );			// semi major axis
	o.P = o.a/fabs(o.a)*2.*M_PI*sqrt(fabs(o.a*o.a*o.a/mu));		// period (negative if hyperbolic)
	o.inc = acos( h2/o.h ) ;				// inclination (wrt xy-plane)   -  Note if pi/2 < i < pi then the orbit is retrograde
	n0 = -h1;					// vector of nodes lies in xy plane => no z component
	n1 =  h0;		
	n = sqrt( n0*n0 + n1*n1 );
	er = p.x*e0 + p.y*e1 + p.z*e2;
	if (n/(o.r*v) <= MIN_REL_ERROR || o.inc <= MIN_REL_ERROR){			// we are in the xy plane
		o.Omega=0.;
		if (o.e <= MIN_REL_ERROR){              // omega not defined for circular orbit
			o.omega = 0.;
		}
		else{
			if (e1>=0.){
				o.omega=acos(e0/o.e);
			}
			else{
				o.omega = 2.*M_PI-acos(e0/o.e);
			}
		}
	}
	else{
		if (o.e <= MIN_REL_ERROR){
			o.omega = 0.;
		}
		else{
			if (e2>=0.){                        // omega=0 if perictr at asc node
				o.omega=acos(( n0*e0 + n1*e1 )/(n*o.e));
			}
			else{
				o.omega=2.*M_PI-acos(( n0*e0 + n1*e1 )/(n*o.e));
			}
		}

		if (n1>=0.){
			o.Omega = acos(n0/n);
		}
		else{
			o.Omega=2.*M_PI-acos(n0/n); // Omega=longitude of asc node
		}								// taken in xy plane from x axis
	}	
	if (o.e<=MIN_REL_ERROR){            // circular orbit
		o.f=0.;                         // f has no meaning
		o.l=0.;
	}
	else{
		double cosf = er/(o.e*o.r);
		double cosea = (1.-o.r/o.a)/o.e;
		
		if (-1.<=cosf && cosf<=1.){     // failsafe
			o.f = acos(cosf);
		}
		else{
			o.f = M_PI/2.*(1.-cosf);
		}
		
		if (-1.<=cosea && cosea<=1.){
			ea  = acos(cosea);
		}
		else{
			ea = M_PI/2.*(1.-cosea);
		}
		
		if (vr<0.){
			o.f=2.*M_PI-o.f;
			ea =2.*M_PI-ea;
		}
		
		o.l = ea -o.e*sin(ea) + o.omega+ o.Omega;  // mean longitude
	}
	return o;
}

void orbit2p(struct particle *p, double G, struct particle *com, struct orbit o){ 
	double r = o.a*(1-o.e*o.e)/(1 + o.e*cos(o.f));

	// Murray & Dermott Eq 2.122
	p->x  = com->x + r*(cos(o.Omega)*cos(o.omega+o.f) - sin(o.Omega)*sin(o.omega+o.f)*cos(o.inc));
	p->y  = com->y + r*(sin(o.Omega)*cos(o.omega+o.f) + cos(o.Omega)*sin(o.omega+o.f)*cos(o.inc));
	p->z  = com->z + r*sin(o.omega+o.f)*sin(o.inc);

	double n = sqrt(G*(p->m+com->m)/(o.a*o.a*o.a));

	// Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velo.incties in the orbital plane
	p->vx = com->vx + (n*o.a/sqrt(1-o.e*o.e))*((o.e+cos(o.f))*(-cos(o.inc)*cos(o.omega)*sin(o.Omega) - cos(o.Omega)*sin(o.omega)) - sin(o.f)*(cos(o.omega)*cos(o.Omega) - cos(o.inc)*sin(o.omega)*sin(o.Omega)));
	p->vy = com->vy + (n*o.a/sqrt(1-o.e*o.e))*((o.e+cos(o.f))*(cos(o.inc)*cos(o.omega)*cos(o.Omega) - sin(o.omega)*sin(o.Omega)) - sin(o.f)*(cos(o.omega)*sin(o.Omega) + cos(o.inc)*cos(o.Omega)*sin(o.omega)));
	p->vz = com->vz + (n*o.a/sqrt(1-o.e*o.e))*((o.e+cos(o.f))*cos(o.omega)*sin(o.inc) - sin(o.f)*sin(o.inc)*sin(o.omega));

	//printf("%f\n", n*n*o.a*o.a/(1.-o.e*o.e)*(1.+o.e*o.e+2.*o.e*cos(o.f))-(p->vx*p->vx + p->vy*p->vy + p->vz*p->vz)); 

}

void move_to_com(struct particle* particles, int N){
	double m = 0;
	double x = 0;
	double y = 0;
	double z = 0;
	double vx = 0;
	double vy = 0;
	double vz = 0;
	for (int i=0;i<N;i++){
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
	for (int i=0;i<N;i++){
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

void forces(struct particle* particles, double t, double dt, double G, int N, int N_megno){	
	struct particle com = particles[0]; // calculate add. forces w.r.t. center of mass
	for(int i=1;i<N;i++){
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
	move_to_com(particles, N);
}

void forces_el(struct particle* particles, double t, double dt, double G, int N, int N_megno){	
	struct particle com = particles[0];
	for(int i=1;i<N;i++){
		struct particle *p = &(particles[i]);
		struct orbit o = tools_p2orbit(G, particles[i], com);
	
		if(i==2){
			printf("%f\n", o.a);
		}
		if (tau_a[i] != 0.){
			double da = -o.a*dt/tau_a[i];
			o.a += da;
			//printf("%f\t%f\n", da, o.a);
		}
		
		if (tau_e[i] != 0.){
			double de = -o.e*dt/tau_e[i];
			o.e += de;
		}

		if (tau_po[i] != 0.){
			double dpo = 2*M_PI*dt/tau_po[i];
			o.omega += dpo;
		}

		orbit2p(&particles[i], G, &com, o); 
		com = get_com(com, particles[i]);
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

	alpha_over_GM0 = Rc/Rc*podot/(1.-gam/2.);
	diskmass = 3.65557*particles[0].m*podot;
	
	particles[0].m += diskmass;
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
	alpha_over_GM0 = 0.;
	diskmass = 0.;

}

