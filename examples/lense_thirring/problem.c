/**
 * Adding gravitational harmonics (J2, J4) to particles
 * 
 * This example shows how to add a J2 and J4 harmonic to particles.
 * If you have GLUT installed for the visualization, press 'w' and/or 'c' for a clearer view of
 * the whole orbit.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    struct reb_particle star = {0};
    star.m     = 1.;
    reb_add(sim, star);


    struct reb_particle planet = {0};  // add a planet on a circular orbit (with default units where G=1)
    planet.x = 1.;
    planet.vy = 1.;
    reb_add(sim, planet);

    struct rebx_extras* rebx = rebx_attach(sim);  // first initialize rebx
    struct rebx_force* lense  = rebx_load_force(rebx, "lense_thirring"); // add our new force
    rebx_add_force(rebx, lense);

    rebx_set_param_double(rebx, &sim->particles[1].ap, "p_hat_x", 0.01);

    double tmax = 100000.;
    reb_integrate(sim, tmax);
}
