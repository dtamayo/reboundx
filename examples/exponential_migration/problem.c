/**
 * Exponential Migration
 *
 * This example shows how to add exponential migration to a REBOUND simulation.
 * Both the Reboundx force (exponential_migration.c) and this example are based on ``modify_orbits_forces''
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* sim);

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Setup constants
    sim->dt             = 0.012;        // initial timestep.
    sim->heartbeat = heartbeat;

    struct reb_particle p = {0}; 
    p.m     = 1.;   
    reb_add(sim, p); 

    double m = 0.0001;
    double a1 = 24.0;
    double e = 0.01;
    double inc = 0.;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;

    struct reb_particle p1 = reb_tools_orbit_to_particle(sim->G, p, m, a1, e, inc, Omega, omega, f);
    reb_add(sim,p1);
    reb_move_to_com(sim);

    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* em = rebx_load_force(rebx, "exponential_migration");  // add force that adds velocity kicks to give exponential migration
    rebx_add_force(rebx, em);

    // Set the timescales for each particle.  
    double tmax = 4.e4;

    rebx_set_param_double(rebx, &sim->particles[1].ap, "em_tau_a", 2.e3);         // add migration e-folding timescale
    rebx_set_param_double(rebx, &sim->particles[1].ap, "em_aini", 24.0);         // add initial semimajor axis
    rebx_set_param_double(rebx, &sim->particles[1].ap, "em_afin", 30.0);         // add final semimajor axis 

	rebx_set_param_int(rebx, &em->ap, "coordinates", REBX_COORDINATES_PARTICLE);
	rebx_set_param_int(rebx, &sim->particles[0].ap, "primary", 1);

    reb_integrate(sim, tmax);
    rebx_free(rebx);    // Free all the memory allocated by rebx
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    // output a e body
    if(reb_output_check(sim, 1.e2)){
        const struct reb_particle sun = sim->particles[0];
        const struct reb_orbit orbit = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sun); // calculate orbit of particles[1]
        printf("%f\t%f\t%f\t%f\n",sim->t,orbit.a, orbit.e);
    }
}
