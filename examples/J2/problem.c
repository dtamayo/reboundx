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
    struct rebx_extras* rebx = rebx_init(sim);

    struct reb_particle star = {0};
    star.m     = 1.;   
    star.hash  = reb_hash("star");
    reb_add(sim, star);

    double m = 0.;
    double a = 1.; 
    double e = 0.2;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle planet = reb_tools_orbit2d_to_particle(sim->G, star, m, a, e, omega, f);
    planet.hash = reb_hash("planet");
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    struct rebx_effect* params = rebx_add(rebx, "gravitational_harmonics");
    double* J2 = rebx_add_param(&sim->particles[1], "J2", REBX_TYPE_DOUBLE);
    double* J4 = rebx_add_param(&sim->particles[1], "J4", REBX_TYPE_DOUBLE);
    double* R_eq = rebx_add_param(&sim->particles[1], "R_eq", REBX_TYPE_DOUBLE);
    *J2 = 10.;
    *J4 = 10.;
    *R_eq = 1/100.;
    
    double tmax = 100.;
    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}
