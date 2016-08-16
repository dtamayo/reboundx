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
    sim->dt = 1.e-8;

    struct reb_particle star = {0};
    star.m     = 1.;   
    star.hash  = reb_hash("star");
    reb_add(sim, star);

    double m = 0.;
    double a = 1.e-4; // put planet close to enhance precession so it's visible in visualization (this would put planet inside the Sun!)
    double e = 0.2;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle planet = reb_tools_orbit2d_to_particle(sim->G, star, m, a, e, omega, f);
    planet.hash = reb_hash("planet");
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    // Could also add "gr" or "gr_full" here.  See documentation for details.
    struct rebx_effect* gr_params = rebx_add(rebx, "gr");
   
    // Have to set speed of light in right units (set by G & initial conditions).  Here we use default units of AU/(yr/2pi)
    double* c = rebx_add_param(gr_params, "c", REBX_TYPE_DOUBLE);  
    *c = REBX_C;

    /* By default, the 'gr' and 'gr_potential' effects assume the massive particle is at index 0.  
     * If it has a different index, or you think the particle might move in the array (e.g. with a custom merger routine)
     * you can add a gr_source flag to the massive particle (for "gr_full" all particles act as sources so this is not necessary):
     */

    int* source = rebx_add_param(&sim->particles[0], "gr_source", REBX_TYPE_INT);
    *source = 1;

    double tmax = 5.e-2;
    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}
