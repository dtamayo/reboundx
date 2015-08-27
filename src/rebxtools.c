#include <math.h>
#include <stdlib.h>
#include "rebound.h"
#include "rebxtools.h"

void rebxtools_orbit2p(double G, struct reb_particle* p, struct reb_particle* primary, struct reb_orbit o){
	int* err = malloc(sizeof(int));
	struct reb_particle p2 = reb_tools_orbit_to_particle_err(G,*primary, p->m, o.a, o.e, o.inc, o.Omega, o.omega, o.f, err);
	p->x = p2.x;
	p->y = p2.y;
	p->z = p2.z;
	p->vx = p2.vx;
	p->vy = p2.vy;
	p->vz = p2.vz;
}
