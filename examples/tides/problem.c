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

int main(int argc, char* argv[]){
	struct reb_simulation* sim = reb_create_simulation();
	// Setup constants
	sim->dt 			= 0.012;		// initial timestep.

	// set up particular parameters
	struct reb_particle p = {0}; 
	p.m  	= 1.;	
	reb_add(sim, p); 

	double m = 0.;
	double a = 1.;
	double e = 0.9;
	double omega = 0.;
	double f = 0.;

	struct reb_particle p2 = reb_tools_orbit2d_to_particle(sim->G, p, m, a, e, omega, f);
	reb_add(sim,p2);

	a = 0.05;
	struct reb_particle p3 = reb_tools_orbit2d_to_particle(sim->G, p, m, a, e, omega, f);
	reb_add(sim,p3);

	reb_move_to_com(sim);

	struct rebx_extras* rebx = rebx_init(sim);	// first initialize rebx

	int N_tides_active = 1; // This makes it so both planets interact tidally with the Sun, but not with each other.  In particular, the first N_tides_active particles will tidally interact with all other particles.  The remaining ones will only interact with the active tide particles (and not each other).  Analogous to sim->N_active for test particles
	rebx_add_tides(rebx, N_tides_active);

	// How to set the various parameters
	
	double Omegax = 0.1; // rotation vector Omega
	double Omegay = 0.2;
	double Omegaz = 1.;
	rebx_set_rot_Omega(&sim->particles[1], Omegax, Omegay, Omegaz);

	double A = 0.8; // principal moments of inertia
	double B = 0.8;
	double C = 1.;
	rebx_set_mom_of_inertia(&sim->particles[1], A, B, C);

	rebx_set_tidal_tau(&sim->particles[1], 0.6); // set tidal time lag
	rebx_set_tidal_k2(&sim->particles[1], 0.01); // set tidal k2

	// mass and radius of objects are set directly
	
	sim->particles[1].m = 1.e-6;
	sim->particles[1].r = 1.e-4;

	// Getting the particle's parameters (you'd use this in tides.c also)
	//
	struct rebx_vec3d Omega = rebx_get_rot_Omega(&sim->particles[1]);
	printf("Particle 1's rotational velocity = (%f, %f, %f)\n", Omega.x, Omega.y, Omega.z);
	struct rebx_mom_of_inertia I = rebx_get_mom_of_inertia(&sim->particles[1]);
	printf("Particle 1's moment of inertia= (%f, %f, %f)\n", I.A, I.B, I.C);

	printf("Particle 1's tidal tau = %f\n", rebx_get_tidal_tau(&sim->particles[1]));
	printf("Particle 1's tidal k2 = %f\n", rebx_get_tidal_k2(&sim->particles[1]));
	printf("Particle 1's mass and radius = (%f, %f)\n", sim->particles[1].m, sim->particles[1].r);

	// set up tide parameters and integrate
	rebx_free(rebx);
}
