/**
 * Galactic tide Example
 *
 * This example shows how to add a galactic tide to a REBOUND simulation.
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

    struct reb_simulation* sim = reb_simulation_create();

    sim->G = 4*M_PI*M_PI;
    reb_simulation_set_integrator(sim, "whfast");
    sim->dt = 10.; 

    sim->heartbeat = heartbeat;

    struct reb_particle star = {0};
    star.m = 1.;
    reb_simulation_add(sim, star);

    double m = 0.;
    double a = 1e5;
    double e = 0;
    double inc = 0;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;

    struct reb_particle particle = reb_particle_from_orbit(sim->G, star, m, a, e, inc, Omega, omega, f);
    reb_simulation_add(sim,particle);
    reb_simulation_move_to_com(sim);

    // add in reboundx and load/add new force
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* gt = rebx_load_force(rebx, "galactic_tide");
    rebx_add_force(rebx, gt);

    double tmax = 1.e8;

    // Set parameters for simulation
    rebx_set_param_double(rebx, &gt->ap, "gt_omega", 26.0); // set galactic angular velocity
    rebx_set_param_double(rebx, &gt->ap, "gt_rho", 0.1);    // set local galactic density

    reb_simulation_integrate(sim, tmax);
    rebx_free(rebx);
    reb_simulation_free(sim);
}

void heartbeat(struct reb_simulation* sim){
    if(reb_simulation_output_check(sim, 1.e6)){
        const struct reb_particle star = sim->particles[0];
        const struct reb_orbit orbit = reb_orbit_from_particle(sim->G, sim->particles[1], star);
        printf("%f\t%f\t%f\t%f\n",sim->t, orbit.a, orbit.e, orbit.inc);
    }
}
