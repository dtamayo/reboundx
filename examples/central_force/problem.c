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

    struct reb_particle star = {0};
    star.m     = 1.;   
    reb_add(sim, star);

    double m = 0.;
    double a = 1.; // put planet close to enhance precession so it's visible in visualization (this would put planet inside the Sun!)
    double e = 0.2;
    double inc = 0.;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle planet = reb_tools_orbit_to_particle(sim->G, star, m, a, e, inc, Omega, omega, f);
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* force = rebx_load_force(rebx, "central_force");
    rebx_add_force(rebx, force);

    /* We first choose a power for our central force (here F goes as r^-1).
     * We then need to add it to the particle(s) that will act as central sources for this force.*/

    struct reb_particle* ps = sim->particles;
    double gammacentral = -1.;
    rebx_set_param_double(rebx, &ps[0].ap, "gammacentral", gammacentral);

    // The other parameter to set is the normalization Acentral (F=Acentral*r^gammacentral). E.g.,

    rebx_set_param_double(rebx, &ps[0].ap, "Acentral", 1.e-4);

    /* We can also use the function rebx_central_force_Acentral to calculate the Acentral required
     * for particles[1] (around primary particles[0]) to have a pericenter precession rate of
     * pomegadot, given a gammacentral value: */
    
    double pomegadot = 1.e-3;
    double Acentral = rebx_central_force_Acentral(ps[1], ps[0], pomegadot, gammacentral);
    rebx_set_param_double(rebx, &ps[0].ap, "Acentral", Acentral);
    
    double tmax = 3.e4;
    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}
