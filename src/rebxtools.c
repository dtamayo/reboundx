#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "reboundx.h"

/* only accepts one reference particle if coordinates=REBX_COORDINATES_PARTICLE.
 * calculate_effect function should check for edge case where particle and reference are the same
 * (could happen e.g. with barycentric coordinates with test particles and single massive body)
 */

struct reb_particle rebx_get_com_without_particle(struct reb_particle com, struct reb_particle p){
    com.x = com.x*com.m - p.x*p.m;
    com.y = com.y*com.m - p.y*p.m;
    com.z = com.z*com.m - p.z*p.m;
    com.vx = com.vx*com.m - p.vx*p.m;
    com.vy = com.vy*com.m - p.vy*p.m;
    com.vz = com.vz*com.m - p.vz*p.m;
    com.ax = com.ax*com.m - p.ax*p.m;
    com.ay = com.ay*com.m - p.ay*p.m;
    com.az = com.az*com.m - p.az*p.m;
    com.m -= p.m;

    if (com.m > 0.){
        com.x /= com.m;
        com.y /= com.m;
        com.z /= com.m;
        com.vx /= com.m;
        com.vy /= com.m;
        com.vz /= com.m;
        com.ax /= com.m;
        com.ay /= com.m;
        com.az /= com.m;
    }
    return com;
}


static inline struct reb_particle rebx_particle_minus(struct reb_particle p1, struct reb_particle p2){
    struct reb_particle p = {0};
    p.m = p1.m-p2.m;
    p.x = p1.x-p2.x;
    p.y = p1.y-p2.y;
    p.z = p1.z-p2.z;
    p.vx = p1.vx-p2.vx;
    p.vy = p1.vy-p2.vy;
    p.vz = p1.vz-p2.vz;
    p.ax = p1.ax-p2.ax;
    p.ay = p1.ay-p2.ay;
    p.az = p1.az-p2.az;
    return p;
}

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

void rebx_com_force(struct reb_simulation* const sim, struct rebx_force* const force, const enum REBX_COORDINATES coordinates, const int back_reactions_inclusive, const char* reference_name, struct reb_vec3d (*calculate_force) (struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p, struct reb_particle* source), struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    struct reb_particle com = reb_get_com(sim); // Start with full com for jacobi and barycentric coordinates.

    int refindex = -1;
    if(coordinates == REBX_COORDINATES_JACOBI){
        refindex = 0;                           // There is no jacobi coordinate for the 0th particle, so set refindex to skip it in loop below.
    }
    else if(coordinates == REBX_COORDINATES_PARTICLE){
        for (int i=0; i < N; i++){
			struct reb_particle* p = &particles[i];
            const int* const reference = rebx_get_param(rebx, p->ap, reference_name);
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
            com = rebx_get_com_without_particle(com, *p);
        }

        struct reb_vec3d a = calculate_force(sim, force, p, &com);
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

void rebxtools_com_ptm(struct reb_simulation* const sim, struct rebx_operator* const operator, const enum REBX_COORDINATES coordinates, const int back_reactions_inclusive, const char* reference_name, struct reb_particle (*calculate_step) (struct reb_simulation* const sim, struct rebx_operator* const operator, struct reb_particle* p, struct reb_particle* source, const double dt), const double dt){
    struct rebx_extras* const rebx = sim->extras;
    const int N_real = sim->N - sim->N_var;
    struct reb_particle com = reb_get_com(sim); // Start with full com for jacobi and barycentric coordinates.

    int refindex = -1;
    if(coordinates == REBX_COORDINATES_JACOBI){
        refindex = 0;                           // There is no jacobi coordinate for the 0th particle, so should skip index 0
    }
    else if(coordinates == REBX_COORDINATES_PARTICLE){
        for (int i=0; i < N_real; i++){
            struct reb_particle* p = &sim->particles[i];
            const int* const reference = rebx_get_param(rebx, p->ap, reference_name);
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
            com = rebx_get_com_without_particle(com, *p);
        }

        struct reb_particle modified_particle = calculate_step(sim, operator, p, &com, dt);
        struct reb_particle diff = rebx_particle_minus(modified_particle, *p);
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

void rebx_tools_calc_plane_Omega_inc(struct reb_vec3d normal_vec, double* Omega, double* inc){
    double r = sqrt(normal_vec.x*normal_vec.x + normal_vec.y*normal_vec.y + normal_vec.z*normal_vec.z);
    *inc = acos(normal_vec.z/r);
    double nx = -normal_vec.y;  // n is the vector pointing toward ascending node
    double ny = normal_vec.x;
    double cosine = nx/sqrt(nx*nx + ny*ny);
    if (cosine > -1. && cosine < 1.){
       *Omega = acos(cosine);
       if (ny < 0){             // if ny is < 0, then azimuthal angle of n is in [-pi, 0]
           *Omega = -*Omega;
       }
    }
    else{
        *Omega = (cosine <= -1.) ? M_PI : 0.;
    }
}

struct reb_vec3d rebx_tools_rotate_XYZ_to_orbital_xyz(struct reb_vec3d XYZ, const double Omega, const double inc, const double omega){
    // adapted from celmech nbody_simulation_utilities.py. Eq 2.121 of 2.8 in Murray & Dermott, xyz = P1^-1P2^-1P3^-1(XYZ)
    double X = XYZ.x;
    double Y = XYZ.y;
    double Z = XYZ.z;

    double s3 = sin(-Omega);
    double c3 = cos(-Omega);
    double x3 = c3 * X - s3 * Y;
    double y3 = s3 * X + c3 * Y;
    double z3 = Z;

    double s2 = sin(-inc);
    double c2 = cos(-inc);
    double x2 = x3;
    double y2 = c2 * y3 - s2 * z3;
    double z2 = s2 * y3 + c2 * z3;

    double s1 = sin(-omega);
    double c1 = cos(-omega);
    double x1 = c1 * x2 - s1 * y2;
    double y1 = s1 * x2 + c1 * y2;
    double z1 = z2;

    struct reb_vec3d xyz = {x1, y1, z1};
    return xyz;
}

struct reb_vec3d rebx_tools_rotate_orbital_xyz_to_XYZ(struct reb_vec3d xyz, const double Omega, const double inc, const double omega){
    // adapted from celmech nbody_simulation_utilities.py. Uses Eq 2.121 of 2.8 in Murray & Dermott, XYZ = P3P2P1(xyz).
    double x = xyz.x;
    double y = xyz.y;
    double z = xyz.z;

    double s1 = sin(omega);
    double c1 = cos(omega);
    double x1 = c1 * x - s1 * y;
    double y1 = s1 * x + c1 * y;
    double z1 = z;

    double s2 = sin(inc);
    double c2 = cos(inc);
    double x2 = x1;
    double y2 = c2 * y1 - s2 * z1;
    double z2 = s2 * y1 + c2 * z1;

    double s3 = sin(Omega);
    double c3 = cos(Omega);
    double x3 = c3 * x2 - s3 * y2;
    double y3 = s3 * x2 + c3 * y2;
    double z3 = z2;

    struct reb_vec3d XYZ = {x3, y3, z3};
    return XYZ;
}

struct reb_vec3d rebx_tools_total_angular_momentum(struct rebx_extras* const rebx){
  struct reb_simulation* const sim = rebx->sim;
	const int N = sim->N;
	const struct reb_particle* restrict const particles = sim->particles;
	const int N_var = sim->N_var;
    struct reb_vec3d L = reb_tools_angular_momentum(sim);
    // Add spin angular momentum of any particles with spin parameters set
    for (int i=0;i<N-N_var;i++){
		struct reb_particle pi = particles[i];
        const double* sx = rebx_get_param(rebx, pi.ap, "sx");
        const double* sy = rebx_get_param(rebx, pi.ap, "sy");
        const double* sz = rebx_get_param(rebx, pi.ap, "sz");
        const double* moi = rebx_get_param(rebx, pi.ap, "moi");

        if (sx != NULL && sy != NULL && sz != NULL && moi != NULL){
          L.x += (*moi) * (*sx);
          L.y += (*moi) * (*sy);
          L.z += (*moi) * (*sz);
        }

	}
	return L;
}

void rebx_align_simulation(struct rebx_extras* const rebx){
    // Modified from celmech nbody_simulation_utilities.py to include spin angular momentum
    struct reb_simulation* const sim = rebx->sim;
    const int N_real = sim->N - sim->N_var;
    double Omega, inc;
    struct reb_vec3d L = rebx_tools_total_angular_momentum(rebx);
    rebx_tools_calc_plane_Omega_inc(L, &Omega, &inc);
    for (int i = 0; i < N_real; i++){
        struct reb_particle* p = &sim->particles[i];
	    struct reb_vec3d pos = {p->x, p->y, p->z};
	    struct reb_vec3d vel = {p->vx, p->vy, p->vz};
        struct reb_vec3d ps = rebx_tools_rotate_XYZ_to_orbital_xyz(pos, Omega, inc, 0);
      	struct reb_vec3d vs = rebx_tools_rotate_XYZ_to_orbital_xyz(vel, Omega, inc, 0);

        // Try this for the spin
        const double* sx = rebx_get_param(rebx, p->ap, "sx");
        const double* sy = rebx_get_param(rebx, p->ap, "sy");
        const double* sz = rebx_get_param(rebx, p->ap, "sz");
        if (sx != NULL && sy != NULL && sz != NULL){
            struct reb_vec3d spin = {*sx, *sy, *sz};
            struct reb_vec3d ss = rebx_tools_rotate_XYZ_to_orbital_xyz(spin, Omega, inc, 0);

            rebx_set_param_double(rebx, &p->ap, "sx", ss.x);
            rebx_set_param_double(rebx, &p->ap, "sy", ss.y);
            rebx_set_param_double(rebx, &p->ap, "sz", ss.z);
        }

      	p->x = ps.x;
      	p->y = ps.y;
      	p->z = ps.z;

      	p->vx = vs.x;
      	p->vy = vs.y;
      	p->vz = vs.z;
  }
}

// TLu transformation matrix
struct reb_vec3d rebx_transform_inv_to_planet(struct reb_vec3d spin_inv, double inc, double omega){
    // This ts a vector from the INVARIANT frame to the PLANET frame
    double sx = spin_inv.x;
    double sy = spin_inv.y;
    double sz = spin_inv.z;

    printf("%e %e %e %e %e\n", sx, sy, sz, inc, omega);

    double t[3][3];

    t[0][0] = cos(omega);
    t[0][1] = sin(omega);
    t[0][2] = 0;
    t[1][0] = -cos(inc) * sin(omega);
    t[1][1] = cos(inc) * cos(omega);
    t[1][2] = sin(inc);
    t[2][0] = sin(inc) * sin(omega);
    t[2][1] = -sin(inc) * cos(omega);
    t[2][2] = cos(inc);

    struct reb_vec3d spin_planet = {0};

    spin_planet.x = sx * t[0][0] + sy * t[0][1] + sz * t[0][2];
    spin_planet.y = sx * t[1][0] + sy * t[1][1] + sz * t[1][2];
    spin_planet.z = sx * t[2][0] + sy * t[2][1] + sz * t[2][2];

    return spin_planet;
}
