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

    rebx_set_param_double(rebx, &sim->particles[0].ap, "lt_rot_rate", 90.361036076);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "lt_Mom_I_fac", 0.06884);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "lt_R_eq", 0.00465247264);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "lt_p_hatx", 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "lt_p_haty", 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "lt_p_hatz", 1.0);

    // Have to set speed of light in right units (set by G & initial conditions).  Here we use default units of AU/(yr/2pi)
    rebx_set_param_double(rebx, &lense->ap, "lt_c", 10065.32);

//  rebx_set_param_double(rebx, &sim->particles[1].ap, "ye_body_density", density);
//  rebx_set_param_int(rebx, &sim->particles[1].ap, "ye_flag", 0);
//  rebx_set_param_double(rebx, &yark->ap, "ye_lstar", lstar);
//  rebx_set_param_double(rebx, &yark->ap, "ye_c", c);


    double tmax = 100000.;
    reb_integrate(sim, tmax);
}
