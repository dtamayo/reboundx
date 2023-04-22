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
    double omega = 90.361036076; //solar rotation rate in rad/year
    double C_I = 0.06884; //solar moment of inertia prefactor
    double R_eq = 0.00465247264; //solar equatorial radius in AU
    double p_x = 0.0; //x-comp. of unit spin-pole direction
    double p_y = 0.0; //y-comp. of unit spin-pole direction
    double p_z = 0.0; //z-comp. of unit spin-pole direction


    struct reb_particle planet = {0};  // add a planet on a circular orbit (with default units where G=1)
    planet.x = 1.;
    planet.vy = 1.;
    reb_add(sim, planet);

    struct rebx_extras* rebx = rebx_attach(sim);  // first initialize rebx
    struct rebx_force* lense  = rebx_load_force(rebx, "lense_thirring"); // add our new force
    rebx_add_force(rebx, lense);

    rebx_set_param_double(rebx, &sim->particles[0].ap, "lt_rot_rate", omega);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "lt_Mom_I_fac", C_I);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "lt_R_eq", R_eq);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "lt_p_hatx", p_x);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "lt_p_haty", p_y);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "lt_p_hatz", p_z);

    // Have to set speed of light in right units (set by G & initial conditions).  Here we use default units of AU/(yr/2pi)
    rebx_set_param_double(rebx, &lense->ap, "lt_c", 10065.32);



    double tmax = 100000.;
    reb_integrate(sim, tmax);
}
