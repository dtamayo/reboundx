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

    //sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->dt = 1.e-8;

    struct reb_particle p = {0}; 
    p.m     = 1.;   
    reb_add(sim, p); 

    double m = 0.;
    double a = 1.e-4; // put planet close to enhance precession so it's visible in visualization (this would put planet inside the Sun!)
    double e = 0.2;
    double omega = 0.;
    double f = 0.;

    struct reb_particle p1 = reb_tools_orbit2d_to_particle(sim->G, p, m, a, e, omega, f);
    reb_add(sim,p1);
    reb_move_to_com(sim);
    
    struct rebx_extras* rebx = rebx_init(sim);
    double c = C_DEFAULT;   // Have to set the speed of light in appropriate units (set by G and your initial conditions).  Here we use the value in default units of AU/(yr/2pi) 
    int source_index = 0;   // Index of the massive particle that is the source of the post-newtonian corrections
    rebx_add_gr(rebx, source_index, c); 
    /*See reboundx.readthedocs.org for more options.*/

    double tmax = 5.e-2;
    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
}
