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
    star.hash  = reb_hash("star");
    reb_add(sim, star);

    double m = 0.;
    double a = 1.; 
    double e = 0.2;
    double inc = 0.;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle planet = reb_tools_orbit_to_particle(sim->G, star, m, a, e, inc, Omega, omega, f);
    planet.hash = reb_hash("planet");
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* gh = rebx_load_force(rebx, "gravitational_harmonics");
    rebx_add_force(rebx, gh);

    rebx_set_param_double(rebx, &sim->particles[0].ap, "J2", 0.1);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "J4", 0.01);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "R_eq", 0.01);
   
    double tmax = 1.e5;
    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}
