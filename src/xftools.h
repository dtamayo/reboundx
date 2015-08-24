#ifndef _XFxftools_H
#define _XFxftools_H

#include "rebound.h"

#ifndef M_PI
#define M_PI 3.1415926535879323846 
#endif

void xftools_inertial_to_jacobi_posvel(void);

struct reb_orbit xftools_orbit_nan(void);

struct reb_orbit xftools_p2orbit(double G, struct reb_particle p, struct reb_particle primary);

void xftools_orbit2p(struct reb_particle* p, double G, struct reb_particle* com, struct reb_orbit o);

void xftools_move_to_com(struct reb_simulation* const r);

struct reb_particle xftools_get_com_of_pair(struct reb_particle p1, struct reb_particle p2);

struct reb_particle xftools_get_com(struct reb_simulation* const r);
#endif
