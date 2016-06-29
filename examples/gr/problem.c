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
    
    sim->particles[0].hash = reb_tools_hash("star");
    sim->particles[1].hash = reb_tools_hash("planet");

    struct rebx_extras* rebx = rebx_init(sim);

    // Could also add "gr" or "gr_full" here.  See documentation for details.
    struct rebx_effect* gr = rebx_add_effect(rebx, "gr_potential");
   
    // Have to set speed of light in right units (set by G & initial conditions).  Here we use default units of AU/(yr/2pi)
    rebx_set_param_double(gr, "c", C_DEFAULT);  
    // Need to set gr_source param to 1 on the massive particle for any effect (except for "gr_full" where all particles act as sources).
    rebx_set_param_int(reb_get_particle_by_name(sim, "star"), "gr_source", 1);
    
    double tmax = 5.e-2;
    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
}
