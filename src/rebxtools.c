#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "rebxtools.h"

static const struct reb_orbit reb_orbit_nan = {.d = NAN, .v = NAN, .h = NAN, .P = NAN, .n = NAN, .a = NAN, .e = NAN, .inc = NAN, .Omega = NAN, .omega = NAN, .pomega = NAN, .f = NAN, .M = NAN, .l = NAN};

#define MIN_REL_ERROR 1.0e-12   ///< Close to smallest relative floating point number, used for orbit calculation
#define TINY 1.E-308        ///< Close to smallest representable floating point number, used for orbit calculation
#define MIN_INC 1.e-8       ///< Below this inclination, the broken angles pomega and theta equal the corresponding 
                            ///< unbroken angles to within machine precision, so a practical boundary for planar orbits
                            //
// returns acos(num/denom), using disambiguator to tell which quadrant to return.  
// will return 0 or pi appropriately if num is larger than denom by machine precision
// and will return 0 if denom is exactly 0.

static double acos2(double num, double denom, double disambiguator){
    double val;
    double cosine = num/denom;
    if(cosine > -1. && cosine < 1.){
        val = acos(cosine);
        if(disambiguator < 0.){
            val = - val;
        }
    }
    else{
        val = (cosine <= -1.) ? M_PI : 0.;
    }
    return val;
}

struct reb_orbit rebxtools_particle_to_orbit_err(double G, struct reb_particle* p, struct reb_particle* primary, int* err){
    *err = 0;
    struct reb_orbit o;
    if (primary->m <= TINY){
        *err = 1;           // primary has no mass.
        return reb_orbit_nan;
    }
    double mu,dx,dy,dz,dvx,dvy,dvz,vsquared,vcircsquared,vdiffsquared;
    double hx,hy,hz,vr,rvr,muinv,ex,ey,ez,nx,ny,n,ea;
    mu = G*(p->m+primary->m);
    dx = p->x - primary->x;
    dy = p->y - primary->y;
    dz = p->z - primary->z;
    dvx = p->vx - primary->vx;
    dvy = p->vy - primary->vy;
    dvz = p->vz - primary->vz;
    o.d = sqrt ( dx*dx + dy*dy + dz*dz );
    
    vsquared = dvx*dvx + dvy*dvy + dvz*dvz;
    o.v = sqrt(vsquared);
    vcircsquared = mu/o.d;  
    o.a = -mu/( vsquared - 2.*vcircsquared );   // semi major axis
    
    hx = (dy*dvz - dz*dvy);                     //angular momentum vector
    hy = (dz*dvx - dx*dvz);
    hz = (dx*dvy - dy*dvx);
    o.h = sqrt ( hx*hx + hy*hy + hz*hz );       // abs value of angular momentum

    vdiffsquared = vsquared - vcircsquared; 
    if(o.d <= TINY){                            
        *err = 2;                                   // particle is on top of primary
        return reb_orbit_nan;
    }
    vr = (dx*dvx + dy*dvy + dz*dvz)/o.d;    
    rvr = o.d*vr;
    muinv = 1./mu;

    ex = muinv*( vdiffsquared*dx - rvr*dvx );
    ey = muinv*( vdiffsquared*dy - rvr*dvy );
    ez = muinv*( vdiffsquared*dz - rvr*dvz );
    o.e = sqrt( ex*ex + ey*ey + ez*ez );        // eccentricity
    o.n = o.a/fabs(o.a)*sqrt(fabs(mu/(o.a*o.a*o.a)));   // mean motion (negative if hyperbolic)
    o.P = 2*M_PI/o.n;                                   // period (negative if hyperbolic)

    o.inc = acos2(hz, o.h, 1.);         // cosi = dot product of h and z unit vectors.  Always in [0,pi], so pass dummy disambiguator
                                        // will = 0 if h is 0.

    nx = -hy;                           // vector pointing along the ascending node = zhat cross h
    ny =  hx;       
    n = sqrt( nx*nx + ny*ny );

    // Omega, pomega and theta are measured from x axis, so we can always use y component to disambiguate if in the range [0,pi] or [pi,2pi]
    o.Omega = acos2(nx, n, ny);         // cos Omega is dot product of x and n unit vectors. Will = 0 if i=0.
    
    ea = acos2(1.-o.d/o.a, o.e, vr);    // from definition of eccentric anomaly.  If vr < 0, must be going from apo to peri, so ea = [pi, 2pi] so ea = -acos(cosea)
    o.M = ea - o.e*sin(ea);                     // mean anomaly (Kepler's equation)

    // in the near-planar case, the true longitude is always well defined for the position, and pomega for the pericenter if e!= 0
    // we therefore calculate those and calculate the remaining angles from them
    if(o.inc < MIN_INC || o.inc > M_PI - MIN_INC){  // nearly planar.  Use longitudes rather than angles referenced to node for numerical stability.
        o.pomega = acos2(ex, o.e, ey);      // cos pomega is dot product of x and e unit vectors.  Will = 0 if e=0.
        o.theta = acos2(dx, o.d, dy);           // cos theta is dot product of x and r vectors (true longitude).  Will = 0 if e = 0.
        if(o.inc < M_PI/2.){
            o.omega = o.pomega - o.Omega;
            o.f = o.theta - o.pomega;
            o.l = o.pomega + o.M;
        }
        else{
            o.omega = o.Omega - o.pomega;
            o.f = o.pomega - o.theta;
            o.l = o.pomega - o.M;
        }
    }
    // in the non-planar case, we can't calculate the broken angles from vectors like above.  omega+f is always well defined, and omega if e!=0
    else{
        double wpf = acos2(nx*dx + ny*dy, n*o.d, dz);   // omega plus f.  Both angles measured in orbital plane, and always well defined for i!=0.
        o.omega = acos2(nx*ex + ny*ey, n*o.e, ez);
        if(o.inc < M_PI/2.){
            o.pomega = o.Omega + o.omega;
            o.f = wpf - o.omega;
            o.theta = o.Omega + wpf;
            o.l = o.pomega + o.M;
        }
        else{
            o.pomega = o.Omega - o.omega;
            o.f = wpf - o.omega;
            o.theta = o.Omega - wpf;
            o.l = o.pomega - o.M;
        }
    }

    return o;
}

struct reb_orbit rebxtools_particle_to_orbit(double G, struct reb_particle* p, struct reb_particle* primary){
    int err;
    struct reb_orbit o = rebxtools_particle_to_orbit_err(G, p, primary, &err);
    if(err != 0){
        fprintf(stderr, "Particle to orbit conversion yielded exception.  Code = %d\n", err);
    }
    return o;
}

void rebxtools_orbit2p(double G, struct reb_particle* p, struct reb_particle* primary, struct reb_orbit* o){
    int err;
    rebxtools_orbit_to_particle(G, p, primary, o->a, o->e, o->inc, o->Omega, o->omega, o->f, &err);
    if(err != 0){
        fprintf(stderr, "Orbit to particle conversion yielded exception.  Code = %d\n", err);
    }
}

static const struct reb_particle reb_particle_nan = {.x = NAN, .y = NAN, .z = NAN, .vx = NAN, .vy = NAN, .vz = NAN, .ax = NAN, .ay = NAN, .az = NAN, .m = NAN, .r = NAN, .lastcollision = NAN, .c = NULL, .id = -1, .ap = NULL, .sim = NULL};

void rebxtools_orbit_to_particle(double G, struct reb_particle* p, struct reb_particle* primary, double a, double e, double inc, double Omega, double omega, double f, int* err){
    *err = 0;
    if(e == 1.){
        *err = 1;       // Can't initialize a radial orbit with orbital elements.
        *p = reb_particle_nan;
        return;
    }
    if(e < 0.){
        *err = 2;       // Eccentricity must be greater than or equal to zero.
        *p = reb_particle_nan;
        return;
    }
    if(e > 1.){
        if(a > 0.){
            *err = 3;   // Bound orbit (a > 0) must have e < 1. 
            *p = reb_particle_nan;
            return;
        }
    }
    else{
        if(a < 0.){
            *err =4;    // Unbound orbit (a < 0) must have e > 1.
            *p = reb_particle_nan;
            return;
        }
    }
    if(e*cos(f) < -1.){
        *err = 5;       // Unbound orbit can't have f set beyond the range allowed by the asymptotes set by the parabola.
        *p = reb_particle_nan;
        return;
    }

    double r = a*(1-e*e)/(1 + e*cos(f));
    double v0 = sqrt(G*(p->m+primary->m)/a/(1.-e*e)); // in this form it works for elliptical and hyperbolic orbits

    double cO = cos(Omega);
    double sO = sin(Omega);
    double co = cos(omega);
    double so = sin(omega);
    double cf = cos(f);
    double sf = sin(f);
    double ci = cos(inc);
    double si = sin(inc);
    
    // Murray & Dermott Eq 2.122
    p->x = primary->x + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
    p->y = primary->y + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
    p->z = primary->z + r*(so*cf+co*sf)*si;

    // Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
    p->vx = primary->vx + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
    p->vy = primary->vy + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
    p->vz = primary->vz + v0*((e+cf)*co*si - sf*si*so);
}

void rebxtools_move_to_com(struct reb_simulation* const sim){
    const int N = sim->N;
    struct reb_particle* restrict const particles = sim->particles;
    struct reb_particle com;
    rebxtools_get_com(sim, sim->N-1, &com);
    for(int i=0; i<N; i++){
        particles[i].x -= com.x;
        particles[i].y -= com.y;
        particles[i].z -= com.z;
        particles[i].vx -= com.vx;
        particles[i].vy -= com.vy;
        particles[i].vz -= com.vz;
    }
}

void rebxtools_update_com_with_particle(struct reb_particle* const com, const struct reb_particle* const p){
    com->x   = com->x*com->m + p->x*p->m;
    com->y   = com->y*com->m + p->y*p->m;
    com->z   = com->z*com->m + p->z*p->m;
    com->vx  = com->vx*com->m + p->vx*p->m;
    com->vy  = com->vy*com->m + p->vy*p->m;
    com->vz  = com->vz*com->m + p->vz*p->m;
    com->m  += p->m;
    if (com->m>0.){
        com->x  /= com->m;
        com->y  /= com->m;
        com->z  /= com->m;
        com->vx /= com->m;
        com->vy /= com->m;
        com->vz /= com->m;
    }
}

void rebxtools_update_com_without_particle(struct reb_particle* const com, const struct reb_particle* const p){
    com->x   = com->x*com->m - p->x*p->m;
    com->y   = com->y*com->m - p->y*p->m;
    com->z   = com->z*com->m - p->z*p->m;
    com->vx  = com->vx*com->m - p->vx*p->m;
    com->vy  = com->vy*com->m - p->vy*p->m;
    com->vz  = com->vz*com->m - p->vz*p->m;
    com->m  -= p->m;
    if (com->m>0.){
        com->x  /= com->m;
        com->y  /= com->m;
        com->z  /= com->m;
        com->vx /= com->m;
        com->vy /= com->m;
        com->vz /= com->m;
    }
}

/* Calculate the center of mass for the first N particles */
void rebxtools_get_com(const struct reb_simulation* const sim, const int first_N, struct reb_particle* com){
    struct reb_particle* particles = sim->particles;
    for (int i=0;i<first_N;i++){
        rebxtools_update_com_with_particle(com, &particles[i]);
    }
}
