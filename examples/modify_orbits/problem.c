/**
 * Migration and other orbit modifications
 *
 * This example shows how to add migration, eccentricity damping
 * and pericenter precession to a REBOUND simulation.  If you have
 * GLUT installed for visualization, press 'w' to see the orbits
 * as wires.  You can zoom out by holding shift, holding down the mouse
 * and dragging.  Press 'c' to better see migration/e-damping.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

double tmax = 5.e4;

int main(int argc, char* argv[]){
	struct reb_simulation* sim = reb_create_simulation();
	// Setup constants
	sim->dt 			= 0.012;		// initial timestep.
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

	// There are two options for how to modify orbits.  You would only choose one (comment the other out).  
	// modify_orbits_forces doesn't have precession implmented yet.

	// modify_orbits_direct directly updates particles' orbital elements at the end of each timestep
	rebx_add_modify_orbits_direct(sim);
	rebx->modify_orbits_direct.tau_a[1] = 1e5;	// add semimajor axis damping on inner planet (e-folding timescale)
	rebx->modify_orbits_direct.tau_e[2] = 1e4;	// add eccentricity damping on outer planet (e-folding timescale)
	rebx->modify_orbits_direct.tau_omega[1] = 1.e4;  // add linear precession on inner planet (precession period)	

	// modify_orbits_forces adds in additional forces that orbit-average to give exponential a and e damping
	/*rebx_add_modify_orbits_forces(sim);
	rebx->modify_orbits_forces.tau_a[1] = 1e5;	// add semimajor axis damping on inner planet (e-folding timescale)
	rebx->modify_orbits_forces.tau_e[2] = 1e4;	// add eccentricity damping on outer planet (e-folding timescale)
	*/

	reb_integrate(sim, tmax);
}
