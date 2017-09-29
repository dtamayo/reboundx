#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "rebxtools.h"
#include "reboundx.h"

/* only accepts one reference particle if coordinates=REBX_COORDINATES_PARTICLE.
 * calculate_effect function should check for edge case where particle and reference are the same
 * (could happen e.g. with barycentric coordinates with test particles and single massive body)
 */

void rebx_calculate_jacobi_masses(const struct reb_particle* const ps, double* const m_j, const int N){
    double eta = ps[0].m;
    for (unsigned int i=1;i<N;i++){ // jacobi masses are reduced mass of particle with interior masses
        m_j[i] = ps[i].m*eta;
        eta += ps[i].m;
        m_j[i] /= eta;
    }
    m_j[0] = eta;
}

double rebx_Edot(struct reb_particle* const ps, const int N){
    double Edot = 0.;
    for(int i=0; i<N; i++){
        Edot += ps[i].m*(ps[i].ax*ps[i].vx + ps[i].ay*ps[i].vy + ps[i].az*ps[i].vz);
    }
    return Edot;
}

void rebx_com_force(struct reb_simulation* const sim, struct rebx_effect* const effect, const enum REBX_COORDINATES coordinates, const int back_reactions_inclusive, const char* reference_name, struct reb_vec3d (*calculate_effect) (struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* p, struct reb_particle* source), struct reb_particle* const particles, const int N){
    struct reb_particle com = reb_get_com(sim); // Start with full com for jacobi and barycentric coordinates.
  
    int refindex = -1;
    if(coordinates == REBX_COORDINATES_JACOBI){
        refindex = 0;                           // There is no jacobi coordinate for the 0th particle, so set refindex to skip it in loop below.
    }
    else if(coordinates == REBX_COORDINATES_PARTICLE){
        for (int i=0; i < N; i++){
			struct reb_particle* p = &particles[i];
            const int* const reference = rebx_get_param_check(p, reference_name, REBX_TYPE_INT);
            if (reference){
                com = particles[i];
                refindex = i;
                break;
            }
            if (i == N-1){
                char str[200];
                sprintf(str, "Coordinates set to REBX_COORDINATES_PARTICLE, but %s param was not found in any particle.  Need to set parameter.\n", reference_name);
                reb_error(sim, str);
            }
        }
    }

    
    for(int i=N-1; i>=0; i--){ // Run through backwards so each iteration does not depend on previous ones in Jacobi coordinates.
        if (i==refindex){
            continue;
        }
        struct reb_particle* p = &particles[i];
        if (coordinates == REBX_COORDINATES_JACOBI){
            com = reb_get_com_without_particle(com, *p);
        }
        
        struct reb_vec3d a = calculate_effect(sim, effect, p, &com); 
        p->ax += a.x;
        p->ay += a.y;
        p->az += a.z;

        double massratio;
        switch(coordinates){
            case REBX_COORDINATES_BARYCENTRIC:
                massratio = p->m/com.m;
                for(int j=0; j < N; j++){
                    particles[j].ax -= massratio*a.x;
                    particles[j].ay -= massratio*a.y;
                    particles[j].az -= massratio*a.z;
                }
                break;
            case REBX_COORDINATES_JACOBI:
                if(back_reactions_inclusive){
                    massratio = p->m/(com.m + p->m);
                }
                else{
                    massratio = p->m/com.m;
                }
                for(int j=0; j < i + back_reactions_inclusive; j++){    // stop at j=i if inclusive, at i-1 if not
                    particles[j].ax -= massratio*a.x;
                    particles[j].ay -= massratio*a.y;
                    particles[j].az -= massratio*a.z;
                }
                break;
            case REBX_COORDINATES_PARTICLE:
                if(back_reactions_inclusive){
                    massratio = p->m/(com.m + p->m);
                    p->ax -= massratio*a.x;
                    p->ay -= massratio*a.y;
                    p->az -= massratio*a.z;
                }
                else{
                    massratio = p->m/com.m;
                }
                particles[refindex].ax -= massratio*a.x;
                particles[refindex].ay -= massratio*a.y;
                particles[refindex].az -= massratio*a.z;
                break;
            default:
                reb_error(sim, "Coordinates not supported in REBOUNDx.\n");
        }
    }
}

static inline void rebx_subtract_posvel(struct reb_particle* p, struct reb_particle* diff, const double massratio){
    p->x -= massratio*diff->x;
    p->y -= massratio*diff->y;
    p->z -= massratio*diff->z;
    p->vx -= massratio*diff->vx;
    p->vy -= massratio*diff->vy;
    p->vz -= massratio*diff->vz;
}

void rebxtools_com_ptm(struct reb_simulation* const sim, struct rebx_effect* const effect, const enum REBX_COORDINATES coordinates, const int back_reactions_inclusive, const char* reference_name, struct reb_particle (*calculate_effect) (struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* p, struct reb_particle* source, const double dt), const double dt){
    const int N_real = sim->N - sim->N_var;
    struct reb_particle com = reb_get_com(sim); // Start with full com for jacobi and barycentric coordinates.
  
    int refindex = -1;
    if(coordinates == REBX_COORDINATES_JACOBI){
        refindex = 0;                           // There is no jacobi coordinate for the 0th particle, so should skip index 0
    }
    else if(coordinates == REBX_COORDINATES_PARTICLE){
        for (int i=0; i < N_real; i++){
            struct reb_particle* p = &sim->particles[i];
            const int* const reference = rebx_get_param_check(p, reference_name, REBX_TYPE_INT);
            if (reference){
                com = sim->particles[i];
                refindex = i;
                break;
            }
            if (i == N_real-1){                 
                char str[200];
                sprintf(str, "Coordinates set to REBX_COORDINATES_PARTICLE, but %s param was not found in any particle.  Need to set parameter.\n", reference_name);
                reb_error(sim, str);
            }
        }
    }

    
    for(int i=N_real-1; i>=0; i--){ // Run through backwards so each iteration does not depend on previous ones in Jacobi coordinates.
        if (i==refindex){
            continue;
        }
        struct reb_particle* p = &sim->particles[i];
        if (coordinates == REBX_COORDINATES_JACOBI){
            com = reb_get_com_without_particle(com, *p);
        }
        
        struct reb_particle modified_particle = calculate_effect(sim, effect, p, &com, dt);
        struct reb_particle diff = reb_particle_minus(modified_particle, *p);
        p->x = modified_particle.x;
        p->y = modified_particle.y;
        p->z = modified_particle.z;
        p->vx = modified_particle.vx;
        p->vy = modified_particle.vy;
        p->vz = modified_particle.vz;

        double massratio;
        switch(coordinates){
            case REBX_COORDINATES_BARYCENTRIC:
                massratio = p->m/com.m;
                for(int j=0; j < N_real; j++){
                    rebx_subtract_posvel(&sim->particles[j], &diff, massratio);
                }
                break;
            case REBX_COORDINATES_JACOBI:
                if(back_reactions_inclusive){
                    massratio = p->m/(com.m + p->m);
                }
                else{
                    massratio = p->m/com.m;
                }
                for(int j=0; j < i + back_reactions_inclusive; j++){    // stop at j=i if inclusive, at i-1 if not
                    rebx_subtract_posvel(&sim->particles[j], &diff, massratio);
                }
                break;
            case REBX_COORDINATES_PARTICLE:
                if(back_reactions_inclusive){
                    massratio = p->m/(com.m + p->m);
                    rebx_subtract_posvel(p, &diff, massratio);
                }
                else{
                    massratio = p->m/com.m;
                }
                rebx_subtract_posvel(&sim->particles[refindex], &diff, massratio);
                break;
            default:
                reb_error(sim, "Coordinates not supported in REBOUNDx.\n");
        }
    }
}
/*static const struct reb_orbit reb_orbit_nan = {.d = NAN, .v = NAN, .h = NAN, .P = NAN, .n = NAN, .a = NAN, .e = NAN, .inc = NAN, .Omega = NAN, .omega = NAN, .pomega = NAN, .f = NAN, .M = NAN, .l = NAN};

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
        val = (cosine <= -1.) ? m_pi : 0.;
    }
    return val;
}

struct reb_orbit rebxtools_particle_to_orbit_err(double g, struct reb_particle* p, struct reb_particle* primary, int* err){
    *err = 0;
    struct reb_orbit o;
    if (primary->m <= tiny){
        *err = 1;           // primary has no mass.
        return reb_orbit_nan;
    }
    double mu,dx,dy,dz,dvx,dvy,dvz,vsquared,vcircsquared,vdiffsquared;
    double hx,hy,hz,vr,rvr,muinv,ex,ey,ez,nx,ny,n,ea;
    mu = g*(p->m+primary->m);
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
    if(o.d <= tiny){                            
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
    o.p = 2*m_pi/o.n;                                   // period (negative if hyperbolic)

    o.inc = acos2(hz, o.h, 1.);         // cosi = dot product of h and z unit vectors.  always in [0,pi], so pass dummy disambiguator
                                        // will = 0 if h is 0.

    nx = -hy;                           // vector pointing along the ascending node = zhat cross h
    ny =  hx;       
    n = sqrt( nx*nx + ny*ny );

    // omega, pomega and theta are measured from x axis, so we can always use y component to disambiguate if in the range [0,pi] or [pi,2pi]
    o.omega = acos2(nx, n, ny);         // cos omega is dot product of x and n unit vectors. will = 0 if i=0.
    
    ea = acos2(1.-o.d/o.a, o.e, vr);    // from definition of eccentric anomaly.  if vr < 0, must be going from apo to peri, so ea = [pi, 2pi] so ea = -acos(cosea)
    o.m = ea - o.e*sin(ea);                     // mean anomaly (kepler's equation)

    // in the near-planar case, the true longitude is always well defined for the position, and pomega for the pericenter if e!= 0
    // we therefore calculate those and calculate the remaining angles from them
    if(o.inc < min_inc || o.inc > m_pi - min_inc){  // nearly planar.  use longitudes rather than angles referenced to node for numerical stability.
        o.pomega = acos2(ex, o.e, ey);      // cos pomega is dot product of x and e unit vectors.  will = 0 if e=0.
        o.theta = acos2(dx, o.d, dy);           // cos theta is dot product of x and r vectors (true longitude).  will = 0 if e = 0.
        if(o.inc < m_pi/2.){
            o.omega = o.pomega - o.omega;
            o.f = o.theta - o.pomega;
            o.l = o.pomega + o.m;
        }
        else{
            o.omega = o.omega - o.pomega;
            o.f = o.pomega - o.theta;
            o.l = o.pomega - o.m;
        }
    }
    // in the non-planar case, we can't calculate the broken angles from vectors like above.  omega+f is always well defined, and omega if e!=0
    else{
        double wpf = acos2(nx*dx + ny*dy, n*o.d, dz);   // omega plus f.  both angles measured in orbital plane, and always well defined for i!=0.
        o.omega = acos2(nx*ex + ny*ey, n*o.e, ez);
        if(o.inc < m_pi/2.){
            o.pomega = o.omega + o.omega;
            o.f = wpf - o.omega;
            o.theta = o.omega + wpf;
            o.l = o.pomega + o.m;
        }
        else{
            o.pomega = o.omega - o.omega;
            o.f = wpf - o.omega;
            o.theta = o.omega - wpf;
            o.l = o.pomega - o.m;
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

void rebxtools_orbit_to_particle(double G, struct reb_particle* p, struct reb_particle* primary, double a, double e, double inc, double Omega, double omega, double f, int* err){
    *err = 0;
    if(e == 1.){
        *err = 1;       // Can't initialize a radial orbit with orbital elements.
        *p = reb_particle_nan();
        return;
    }
    if(e < 0.){
        *err = 2;       // Eccentricity must be greater than or equal to zero.
        *p = reb_particle_nan();
        return;
    }
    if(e > 1.){
        if(a > 0.){
            *err = 3;   // Bound orbit (a > 0) must have e < 1. 
            *p = reb_particle_nan();
            return;
        }
    }
    else{
        if(a < 0.){
            *err =4;    // Unbound orbit (a < 0) must have e > 1.
            *p = reb_particle_nan();
            return;
        }
    }
    if(e*cos(f) < -1.){
        *err = 5;       // Unbound orbit can't have f set beyond the range allowed by the asymptotes set by the parabola.
        *p = reb_particle_nan();
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

void rebxtools_get_com(const struct reb_simulation* const sim, const int first_N, struct reb_particle* com){
    struct reb_particle* particles = sim->particles;
    for (int i=0;i<first_N;i++){
        rebxtools_update_com_with_particle(com, &particles[i]);
    }
}*/
