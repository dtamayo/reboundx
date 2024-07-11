/**
 * Gas Damping Forces Example
 *
 * This example shows how to add gas damping forces to a REBOUND simulation.
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

    /* start the rebound simulation here */
    struct reb_simulation* sim = reb_simulation_create();

    sim->heartbeat = heartbeat;

    // add the host star
    struct reb_particle star = {0}; 
    star.m     = 1.;   
    reb_simulation_add(sim, star); 

    double m = 3.e-6;          // roughly 1 Earth Mass in Solar Masses
    double a = 0.1;
    double e = 0.05;
    double inc = 5.*M_PI/180.;        // 5 degrees in radians
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;

    struct reb_particle planet = reb_particle_from_orbit(sim->G, star, m, a, e, inc, Omega, omega, f);
    reb_simulation_add(sim,planet);

    // move simulation to center-of-mass frame
    reb_simulation_move_to_com(sim);

    // add in reboundx and load/add new force
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* gd = rebx_load_force(rebx, "gas_damping_forces");
    rebx_add_force(rebx, gd);

    // Set the total simulation time  
    double tmax = 5.e2*2*M_PI;      // 500 years in yr/2pi

    // Set parameter for particle
    rebx_set_param_double(rebx, &sim->particles[1].ap, "d_factor", 10.);         // set depletion factor to 10

    // integrate simulation
    reb_simulation_integrate(sim, tmax);
    rebx_free(rebx);            // free all the memory allocated by reboundx
    reb_simulation_free(sim);   // free all the memory allocated by rebound 
}

void heartbeat(struct reb_simulation* sim){
    // output a e i of the planet
    if(reb_simulation_output_check(sim, 50.*2*M_PI)){
        const struct reb_particle star = sim->particles[0];
        const struct reb_orbit orbit = reb_orbit_from_particle(sim->G, sim->particles[1], star); // calculate orbit of planet
        printf("%f\t%f\t%f\t%f\n",sim->t, orbit.a, orbit.e, orbit.inc);
    }
}
