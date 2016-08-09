/**
 * Precession from tides.
 *
 * This example shows how to add precession due to tides raised on either the primary, the orbiting bodies, or both.
 * See also the corresponding Python example for more details.
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
    sim->dt = 1.e-5;

    struct reb_particle star = {0};
    star.m     = 1.;   
    star.hash  = reb_hash("star");
    reb_add(sim, star);

    double m = 1.e-5; // if planets don't have mass, they won't feel tides.
    double a = 3.e-2; // put planet close to enhance precession so it's visible in visualization (this would put planet inside the Sun!)
    double e = 0.2;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle planet = reb_tools_orbit2d_to_particle(sim->G, star, m, a, e, omega, f);
    planet.hash = reb_hash("planet");
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    struct rebx_effect* tides = rebx_add(rebx, "tides_precession");
   
    // Have to set particle radii and k2 Love numbers
    sim->particles[0].r = 0.005; // in consistent units (here we're using default G=1, so AU)
    sim->particles[1].r = 0.0005;
    rebx_set_particle_param_double(&sim->particles[0], "k2", 0.03);
    rebx_set_particle_param_double(&sim->particles[1], "k2", 0.3);

    /* By default, implementation assumes particles[0] is the primary.
     * You can also set the primary flag explicitly:
     */

    rebx_set_particle_param_int(&sim->particles[0], "primary", 1);

    double tmax = 1000.;
    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}
