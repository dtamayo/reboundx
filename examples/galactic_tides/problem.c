/**
 * Galactic Tides Example
 *
 * This example shows how to add galactic tides to a REBOUND simulation.
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

    // set units to use throughout the simulation
    sim->G = 4*M_PI*M_PI;  //  Gravitational constant in AU^3/M_sun/yr^2

    sim->heartbeat = heartbeat;

    // add the host star
    struct reb_particle star = {0}; 
    star.m     = 1.;   
    reb_simulation_add(sim, star); 

    double m = 0.;                // test particle
    double a = 1e5;
    double e = 1-25./a;
    double inc = 1.*M_PI/180.;        // 1 degree in radians
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;

    struct reb_particle particle = reb_particle_from_orbit(sim->G, star, m, a, e, inc, Omega, omega, f);
    reb_simulation_add(sim,particle);

    // move simulation to center-of-mass frame
    reb_simulation_move_to_com(sim);

    // add in reboundx and load/add new force
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* gt = rebx_load_force(rebx, "galactic_tides");
    rebx_add_force(rebx, gt);

    // Set the total simulation time  
    double tmax = 4.e9;      // 4 billion years

    // Set parameters for simulation
    rebx_set_param_double(rebx, &gt->ap, "gt_omega", 26.0);              // set galactic angular velocity
    rebx_set_param_double(rebx, &gt->ap, "gt_rho", 0.1);                // set local galactic density

    // integrate simulation
    reb_simulation_integrate(sim, tmax);
    rebx_free(rebx);            // free all the memory allocated by reboundx
    reb_simulation_free(sim);   // free all the memory allocated by rebound 
}

void heartbeat(struct reb_simulation* sim){
    // output a e i of the planet
    if(reb_simulation_output_check(sim, 50.)){
        const struct reb_particle star = sim->particles[0];
        const struct reb_orbit orbit = reb_orbit_from_particle(sim->G, sim->particles[1], star); // calculate orbit of planet
        printf("%f\t%f\t%f\t%f\n",sim->t, orbit.a, orbit.e, orbit.inc);
    }
}
