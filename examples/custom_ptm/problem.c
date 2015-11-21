/**
 * Adding custom post-timestep modifications and forces.
 *
 * This allows the user to use the built-in functions of REBOUNDx
 * but also include their own specialised functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

// A very simple post timestep modification to change the planet's orbit
void simple_drag(struct reb_simulation* const r){
	r->particles[2].vx *= 1. - 1e-5*r->dt_last_done;
	r->particles[2].vy *= 1. - 1e-5*r->dt_last_done;
	r->particles[2].vz *= 1. - 1e-5*r->dt_last_done;
}

int main(int argc, char* argv[]){
	struct reb_simulation* sim = reb_create_simulation();
	// Setup constants
	sim->dt 			= 0.012;	// initial timestep.
	sim->integrator	= REB_INTEGRATOR_IAS15;

	struct reb_particle p = {0}; 
	p.m  	= 1.;	
	reb_add(sim, p); 

	double m = 0.;
	double a1 = 1.;
	double a2 = 2.;
	double e = 0.4;
	double omega = 0.;
	double f = 0.;

	struct reb_particle p1 = reb_tools_orbit2d_to_particle(sim->G, p, m, a1, e, omega, f);
	struct reb_particle p2 = reb_tools_orbit2d_to_particle(sim->G, p, m, a2, e, omega, f);
	reb_add(sim,p1);
	reb_add(sim,p2);
	reb_move_to_com(sim);

	struct rebx_extras* rebx = rebx_init(sim);	// first initialize rebx

	rebx_add_custom_ptm(rebx,simple_drag);		// Add additional post timestep modifications

	double tmax = 5.e4;
	reb_integrate(sim, tmax);
	rebx_free(rebx);				// Free all the memory allocated by rebx
}
