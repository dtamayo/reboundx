/**
 * Post-Newtonian correction from general relativity
 * 
 * This example shows how to add post-newtonian corrections to REBOUND simulations with reboundx.
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
    double e = 0.5;
    double omega = 0.;
    double f = M_PI;
    
    struct reb_particle planet = reb_tools_orbit2d_to_particle(sim->G, star, m, a, e, omega, f);
    planet.hash = reb_hash("planet");
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    struct rebx_effect* track_min_distance = rebx_add(rebx, "track_min_distance");
   
    double* min_distance = rebx_add_param(&sim->particles[1], "min_distance", REBX_TYPE_DOUBLE);
    *min_distance = 5.;

    double tmax = 10.;
    reb_integrate(sim, tmax); 
    min_distance = rebx_get_param(&sim->particles[1], "min_distance");
    printf("%e\n", *min_distance);
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}
