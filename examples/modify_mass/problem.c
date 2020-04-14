/**
 * Exponential mass loss/gain
 *
 * This example shows how to add exponential mass loss to bodies in the integration.
 * The primary's mass loss causes planets' orbits to expand.
 * For more flexible variations in mass or other parameters see the parameter interpolation example.
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
	
	struct rebx_extras* rebx = rebx_attach(sim); // initialize reboundx
    struct rebx_operator* modify_mass = rebx_load_operator(rebx, "modify_mass");

    /* The function rebx_add_operator will choose how to add the operator to the integration
     * scheme based on the integrator being used and the properties of the operator.
     * This is typically a half operator timestep before the main REBOUND timestep, and half afterward.
     */

	rebx_add_operator(rebx, modify_mass);

    /* If you wanted to make your own choices, you can add individual operator steps.
     * In this case you would pass additional parameters. Say we wanted to add a full operaator timestep after the main REBOUND timestep;
     *
     * dt_fraction = 1. // Fraction of a REBOUND timestep (sim->dt) operator should act
     * timing = REBX_TIMING_POST; // Should happen POST timestep
     * name = "modify_mass_post"; // Name identifier
     * rebx_add_operator_step(rebx, modify_mass, dt_fraction, timing, name);
     */

	// To set an exponential mass loss rate, we set the e-folding timescale (positive for growth, negative for loss)
    // Here have the star lose mass with e-damping timescale = tmax
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tau_mass", -tmax);

	// We can approximate a linear mass loss/growth rate if the rate is small by taking tau_mass = M_initial / mass_loss_rate (or growth)
	double M_dot = 1.e-12; // mass growth rate for the planet (in simulation units--here Msun/yr)
    double tau_mass = sim->particles[1].m / M_dot; // first planet gains mass at linear rate M_dot
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_mass", tau_mass);

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
