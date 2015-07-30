#include <math.h>
#include "xftools.h"

/*void xftools_inertial_to_jacobi_posvel(void){
	double s_x = particles[N_megno].m * particles[N_megno].x;
	double s_y = particles[N_megno].m * particles[N_megno].y;
	double s_z = particles[N_megno].m * particles[N_megno].z;
	double s_vx = particles[N_megno].m * particles[N_megno].vx;
	double s_vy = particles[N_megno].m * particles[N_megno].vy;
	double s_vz = particles[N_megno].m * particles[N_megno].vz;
	for (unsigned int i=1+N_megno;i<N;i++){
		const double ei = etai[i-1-N_megno];
		const struct particle pi = particles[i];
		const double pme = eta[i-N_megno]*ei;
		p_j[i].x = pi.x - s_x*ei;
		p_j[i].y = pi.y - s_y*ei;
		p_j[i].z = pi.z - s_z*ei;
		p_j[i].vx = pi.vx - s_vx*ei;
		p_j[i].vy = pi.vy - s_vy*ei;
		p_j[i].vz = pi.vz - s_vz*ei;
		s_x  = s_x  * pme + pi.m*p_j[i].x ;
		s_y  = s_y  * pme + pi.m*p_j[i].y ;
		s_z  = s_z  * pme + pi.m*p_j[i].z ;
		s_vx = s_vx * pme + pi.m*p_j[i].vx;
		s_vy = s_vy * pme + pi.m*p_j[i].vy;
		s_vz = s_vz * pme + pi.m*p_j[i].vz;
	}
	p_j[N_megno].x = s_x * Mtotali;
	p_j[N_megno].y = s_y * Mtotali;
	p_j[N_megno].z = s_z * Mtotali;
	p_j[N_megno].vx = s_vx * Mtotali;
	p_j[N_megno].vy = s_vy * Mtotali;
	p_j[N_megno].vz = s_vz * Mtotali;
}*/

struct reb_orbit xftools_orbit_nan(void){
	struct reb_orbit o;
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

struct reb_orbit xftools_p2orbit(double G, struct reb_particle p, struct reb_particle primary){
	struct reb_orbit o;
	double h0,h1,h2,e0,e1,e2,n0,n1,n,er,vr,mu,ea,dx,dy,dz,dvx,dvy,dvz,v,cosf,cosea;
	mu = G*(p.m+primary.m);
	dx = p.x - primary.x;
	dy = p.y - primary.y;
	dz = p.z - primary.z;
	dvx = p.vx - primary.vx;
	dvy = p.vy - primary.vy;
	dvz = p.vz - primary.vz;
	h0 = (dy*dvz - dz*dvy); 			//angular momentum vector
	h1 = (dz*dvx - dx*dvz);
	h2 = (dx*dvy - dy*dvx);
	o.h = sqrt ( h0*h0 + h1*h1 + h2*h2 );		// abs value of angular moment 
	v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
	o.r = sqrt ( dx*dx + dy*dy + dz*dz );
	if (o.h/(o.r*v) <= MIN_REL_ERROR){
		return xftools_orbit_nan();
	}
	vr = (dx*dvx + dy*dvy + dz*dvz)/o.r;
	e0 = 1./mu*( (v*v-mu/o.r)*dx - o.r*vr*dvx );
	e1 = 1./mu*( (v*v-mu/o.r)*dy - o.r*vr*dvy );
	e2 = 1./mu*( (v*v-mu/o.r)*dz - o.r*vr*dvz );
 	o.e = sqrt( e0*e0 + e1*e1 + e2*e2 );		// eccentricity
	o.a = -mu/( v*v - 2.*mu/o.r );			// semi major axis

	o.P = o.a/fabs(o.a)*2.*M_PI*sqrt(fabs(o.a*o.a*o.a/mu));		// period (negative if hyperbolic)
	o.inc = acos( h2/o.h ) ;				// inclination (wrt xy-plane)   -  Note if pi/2 < i < pi then the orbit is retrograde
	n0 = -h1;					// vector of nodes lies in xy plane => no z component
	n1 =  h0;		
	n = sqrt( n0*n0 + n1*n1 );
	er = dx*e0 + dy*e1 + dz*e2;
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
		cosf = er/(o.e*o.r);
		cosea = (1.-o.r/o.a)/o.e;
		
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

void xftools_orbit2p(struct reb_particle* p, double G, struct reb_particle* com, struct reb_orbit o){ 
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

void xftools_move_to_com(struct reb_particle* particles, int N){
	double m = 0;
	double x = 0;
	double y = 0;
	double z = 0;
	double vx = 0;
	double vy = 0;
	double vz = 0;
	for (int i=0;i<N;i++){
		struct reb_particle p = particles[i];
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

struct reb_particle xftools_get_com(struct reb_particle p1, struct reb_particle p2){
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
