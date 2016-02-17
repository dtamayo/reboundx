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

void heartbeat(struct reb_simulation* sim);
double tmax = 1.e4;

int main(int argc, char* argv[]){
	struct reb_simulation* sim = reb_create_simulation();
    sim->G = 4*M_PI*M_PI;               // use units of AU, yr and solar masses
	sim->heartbeat = heartbeat;
	
	struct reb_particle sun = {0}; 
	sun.m  	= 1.;	
	reb_add(sim, sun); 

	struct reb_particle planet = {0};   // Initialize planets on circular orbits, each 2 times farther than last.
	planet.m = 1.e-3;
	planet.x = 1.;
	planet.vy = 2.*M_PI;
	reb_add(sim,planet);
	planet.x *= 2.;
	planet.vy /= sqrt(2.);
	reb_add(sim, planet);
	planet.x *= 2.;
	planet.vy /= sqrt(2.);
	reb_add(sim, planet);

	reb_move_to_com(sim);
	
	struct rebx_extras* rebx = rebx_init(sim); // initialize reboundx
	rebx_add_modify_mass(rebx); 

	// To set an exponential mass loss rate, we set the e-folding timescale (positive for growth, negative for loss)
    rebx_set_param_double(&sim->particles[0], "tau_mass", -tmax); // star loses mass with e-damping timescale tmax

	// We can approximate a linear mass loss/growth rate if the rate is small by taking tau_mass = M_initial / mass_loss_rate (or growth)
	double M_dot = 1.e-12; // mass growth rate for the planet (in simulation units--here Msun/yr)
    double tau_mass = sim->particles[1].m / M_dot;
    rebx_set_param_double(&sim->particles[1], "tau_mass", tau_mass); // first planet gains mass at linear rate M_dot

	reb_integrate(sim, tmax); 
	rebx_free(rebx); 	// this explicitly frees all the memory allocated by REBOUNDx 
}

void heartbeat(struct reb_simulation* const sim){ 
	// Output masses and semimajor of the inner planet 100 times over the time of the simulation
    if(reb_output_check(sim, tmax/100.)){
		struct reb_orbit o = rebxtools_particle_to_orbit(sim->G, &sim->particles[1], &sim->particles[0]);
		printf("t=%e, Sun mass = %f, planet mass = %e, planet semimajor axis = %f\n", sim->t, sim->particles[0].m, sim->particles[1].m, o.a);
	}   
 }
