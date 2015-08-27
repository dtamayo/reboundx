#ifndef _REBXTOOLS_H
#define _REBXTOOLS_H

#include "rebound.h"

#ifndef M_PI
#define M_PI 3.1415926535879323846 
#endif

void rebxtools_orbit2p(double G, struct reb_particle* p, struct reb_particle* primary, struct reb_orbit o);
#endif
