#include "gr.h"
#include "reboundxf.h"
#include "rebound.h"
#include <stdio.h>

void rebxf_gr(struct reb_simulation* const sim){
	struct rebxf_params* xf = sim->xf_params;
	struct rebxf_param_gr xfparams = xf->gr;
	const double C = xfparams.c;
	const int _N_real = sim->N - sim->N_var;
	const double G = sim->G;
	struct reb_particle* const particles = sim->particles;
	
	for (int i=1; i<_N_real; i++){
		const struct reb_particle sun = particles[0];
		const double dx = particles[i].x - sun.x;
		const double dy = particles[i].y - sun.y;
		const double dz = particles[i].z - sun.z;
		const double r2 = dx*dx + dy*dy + dz*dz;
		const double _r = sqrt(r2);
		const double dvx = particles[i].vx - sun.vx;
		const double dvy = particles[i].vy - sun.vy;
		const double dvz = particles[i].vz - sun.vz;
		// Benitez and Gallardo 2008
		const double alpha = G*sun.m/(_r*_r*_r*C*C);
		const double v_mag = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
		const double beta = 4.*G*sun.m/_r - v_mag*v_mag;
		const double gamma = 4.*(dvx*dx + dvy*dy + dvz*dz);

		const double dax = alpha*(beta*dx + gamma*dvx);
		const double day = alpha*(beta*dy + gamma*dvy);
		const double daz = alpha*(beta*dz + gamma*dvz);

		particles[i].ax += dax;
		particles[i].ay += day;
		particles[i].az += daz;
		particles[0].ax -= particles[i].m/particles[0].m*dax;
		particles[0].ay -= particles[i].m/particles[0].m*day;
		particles[0].az -= particles[i].m/particles[0].m*daz;
	}
}

