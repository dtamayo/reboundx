#ifndef _XFxftools_H
#define _XFxftools_H

#include "libreboundxf.h"

#ifndef M_PI
#define M_PI 
#endif

void xftools_inertial_to_jacobi_posvel(void);

struct orbit xftools_orbit_nan(void);

struct orbit xftools_p2orbit(double G, struct particle p, struct particle primary);

void xftools_orbit2p(struct particle *p, double G, struct particle *com, struct orbit o);

void xftools_move_to_com(struct particle* particles, int N);

struct particle xftools_get_com(struct particle p1, struct particle p2);

#endif
