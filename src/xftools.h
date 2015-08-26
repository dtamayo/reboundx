#ifndef _XFxftools_H
#define _XFxftools_H

#include "rebound.h"

#ifndef M_PI
#define M_PI 3.1415926535879323846 
#endif

void xftools_inertial_to_jacobi_posvel(void);

struct reb_orbit xftools_orbit_nan(void);

struct reb_orbit xftools_particle_to_orbit(double G, struct reb_particle p, struct reb_particle primary, int* err);

struct reb_particle xftools_orbit_to_particle(double G, struct reb_particle primary, double m, double a, double e, double inc, double Omega, double omega, double f, int* err);

void xftools_orbit2p(double G, struct reb_particle* p, struct reb_particle* primary, struct reb_orbit o);

void xftools_move_to_com(struct reb_simulation* const r);

struct reb_particle xftools_get_com_of_pair(struct reb_particle p1, struct reb_particle p2);

struct reb_particle xftools_get_com(struct reb_simulation* const r);
#endif
