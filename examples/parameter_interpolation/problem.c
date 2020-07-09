/**
 * Stellar evolution with interpolated mass data
 *
 * This example shows how to change a particle's mass by interpolating
 * time-series data during a REBOUND simulation. If you have GLUT installed for
 * visualization, press 'w' to see the orbits as wires. You can zoom out by
 * holding shift, holding down the mouse and dragging. Press 'c' to better see
 * migration/e-damping.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* sim);
struct rebx_interpolator* stellarmass;
struct rebx_extras* rebx;
double tmax = 1e4;
int main(int argc, char* argv[]) {
	struct reb_simulation* sim = reb_create_simulation();
    sim->G = 4*M_PI*M_PI; // use units of AU, yr and solar masses
	sim->heartbeat = heartbeat;
	sim->integrator = REB_INTEGRATOR_WHFAST;
	
	struct reb_particle sun = {0};
	sun.m  	= 1.;
	reb_add(sim, sun);
	// Initialize planets on circular orbits, each 2 times farther than last.
	struct reb_particle planet = {0};
	planet.x = 1.;
	planet.vy = 2.*M_PI;
	reb_add(sim, planet);
	planet.x *= 2.;
	planet.vy /= sqrt(2.);
	reb_add(sim, planet);
	planet.x *= 2.;
	planet.vy /= sqrt(2.);
	reb_add(sim, planet);

	reb_move_to_com(sim);
	rebx = rebx_attach(sim); // initialize reboundx

	// To set how the mass will change, we pass three equally-sized arrays
	// necessary to interpolate the time-mass values. Here we have six (6)
	// values that correspond to a star losing mass with an e-damping timescale
	// of tmax (1e4) up to 12,500 yr. The effect will use a cubic spline to
	// interpolate any intermediate values needed by the simulation.
    int n = 6; // size of arrays
	double times[] = {0, 2500, 5000, 7500, 10000, 12500}; // in yr
	double values[] = {1., 0.77880078307, 0.60653065971, 0.47236655274, 0.36787944117, 0.28650479686}; // in Msun
    stellarmass = rebx_create_interpolator(rebx, n, times, values, REBX_INTERPOLATION_SPLINE);

    reb_integrate(sim, tmax); 
    rebx_free_interpolator(stellarmass);
	rebx_free(rebx); // explicitly free all the memory allocated by REBOUNDx
}

void heartbeat(struct reb_simulation* sim) {
	if (reb_output_check(sim, tmax/1000)) {			
        sim->particles[0].m = rebx_interpolate(rebx, stellarmass, sim->t);
        reb_move_to_com(sim);
        struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
        printf("t=%e, Sun mass = %f, planet semimajor axis = %f\n", sim->t, sim->particles[0].m, o.a);
	}
}
