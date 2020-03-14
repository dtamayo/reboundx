/**
 * Stellar evolution with splined mass data
 *
 * This example shows how to change a particle's mass using splined time-series data during a REBOUND simulation. 
 * If you have GLUT installed for visualization, press 'w' to see the orbits as wires.
 * You can zoom out by holding shift, holding down the mouse and dragging.
 * Press 'c' to better see migration/e-damping.
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
    sim->G = 4*M_PI*M_PI;                   // use units of AU, yr and solar masses
	sim->heartbeat = heartbeat;

	sim->integrator = REB_INTEGRATOR_IAS15;
	sim->ri_ias15.epsilon = 0;              // makes IAS15 non-adaptive
	
	struct reb_particle sun = {0}; 
	sun.m  	= 1.;	
	reb_add(sim, sun); 

	struct reb_particle planet = {0};       // Initialize planets on circular orbits, each 2 times farther than last.
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
	
	struct rebx_extras* rebx = rebx_attach(sim); // initialize reboundx
    struct rebx_operator* stellar_evo = rebx_load_operator(rebx, "stellar_evo");

    /* The function rebx_add_operator will choose how to add the operator to the integration
     * scheme based on the integrator being used and the properties of the operator.
     * This is typically a half operator timestep before the main REBOUND timestep, and half afterward.
     */
	rebx_add_operator(rebx, stellar_evo);

    /* If you wanted to make your own choices, you can add individual operator steps.
     * In this case you would pass additional parameters. Say we wanted to add a full operator timestep after the main REBOUND timestep;
     *
     * dt_fraction = 1. // Fraction of a REBOUND timestep (sim->dt) operator should act
     * timing = REBX_TIMING_POST; // Should happen POST timestep
     * name = "stellar_evo_post"; // Name identifier
     * rebx_add_operator_step(rebx, stellar_evo, dt_fraction, timing, name);
     */

	// To set how the mass will change, we pass three equally-sized arrays necessary to interpolate the time-mass values.
    // Here we have six (6) values that correspond to a star losing mass with an e-damping timescale of -tmax (-1e4) up to 12,500 yr.
	// The effect will use a cubic spline to interpolate any intermediate values needed by the simulation.
    int n = 6;																						     // size of arrays
	double mass_age[] = {0, 2500, 5000, 7500, 10000, 12500}; 											 // in yr
	double mass_val[] = {1., 0.77880078307, 0.60653065971, 0.47236655274, 0.36787944117, 0.28650479686}; // in Msun
    double mass_2val[n];                                                                                 // always uninitialized, n-sized array

    rebx_set_param_int(rebx, &sim->particles[0].ap, "mass_n", n);
	rebx_set_param_pointer(rebx, &sim->particles[0].ap, "mass_age", mass_age);
	rebx_set_param_pointer(rebx, &sim->particles[0].ap, "mass_val", mass_val);
    rebx_set_param_pointer(rebx, &sim->particles[0].ap, "mass_2val", mass_2val);
	
	reb_integrate(sim, tmax); 
	rebx_free(rebx); 	// this explicitly frees all the memory allocated by REBOUNDx 
}

void heartbeat(struct reb_simulation* const sim){ 
	// Output masses and semimajor of the inner planet 100 times over the time of the simulation
    if(reb_output_check(sim, tmax/100.)){
		struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
		printf("t=%e, Sun mass = %f, planet mass = %e, planet semimajor axis = %f\n", sim->t, sim->particles[0].m, sim->particles[1].m, o.a);
	}   
 }
