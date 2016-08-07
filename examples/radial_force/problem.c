/**
 * General central force.
 * 
 * This example shows how to add a general central force.
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
    reb_add(sim, star);

    double m = 0.;
    double a = 1.; // put planet close to enhance precession so it's visible in visualization (this would put planet inside the Sun!)
    double e = 0.2;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle planet = reb_tools_orbit2d_to_particle(sim->G, star, m, a, e, omega, f);
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    rebx_add(rebx, "radial_force");
    /* We first choose a power (must be a double!) for our radial force (here F goes as r^-1).
     * We then need to add it to the particle(s) that will act as central sources for this force.*/

    double gammaradial = -1.; // must be a double
    
    struct reb_particle* ps = sim->particles;
    rebx_set_param_double(&ps[0], "gammaradial", gammaradial);

    // The other parameter to set is the normalization Aradial (F=Aradial*r^gammaradial). E.g.,

    rebx_set_param_double(&ps[0], "Aradial", 1.e-4);

    /* We can also use the function rebx_radial_force_Aradial to calculate the Aradial required
     * for particles[1] (around primary particles[0]) to have a pericenter precession rate of
     * pomegadot, given a gammaradial value: */
    
    double pomegadot = 1.e-3;
    double Aradial = rebx_radial_force_Aradial(ps[1], ps[0], pomegadot, gammaradial);
    rebx_set_param_double(&ps[0], "Aradial", Aradial);
    
    double tmax = 3.e4;
    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}
