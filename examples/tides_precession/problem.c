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

    double m = 1.e-5;   // if planets don't have mass, they won't feel tides.
    double a = 3.e-2;   // put planet close to enhance precession so it's visible in visualization (this would put planet inside the Sun!)
    double e = 0.2;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle planet = reb_tools_orbit2d_to_particle(sim->G, star, m, a, e, omega, f);
    planet.hash = reb_hash("planet");
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    struct rebx_effect* tides = rebx_add(rebx, "tides_precession");
   
                        // Have to set R_tides (physical radius) and k1 (apsidal motion constant, half the tidal Love number k2) parameters. 
                        // Could just set one set to consider tides raised on one body only.
    double* R_tides = rebx_add_param(&sim->particles[0], "R_tides", REBX_TYPE_DOUBLE);
    *R_tides = 0.005;         // in consistent units (here we're using default G=1, so AU)
    R_tides = rebx_add_param(&sim->particles[1], "R_tides", REBX_TYPE_DOUBLE);
    *R_tides = 0.0005;
    double* k1 = rebx_add_param(&sim->particles[0], "k1", REBX_TYPE_DOUBLE);
    *k1 = 0.03;
    k1 = rebx_add_param(&sim->particles[1], "k1", REBX_TYPE_DOUBLE);
    *k1 = 0.3;

    /* By default, implementation assumes particles[0] is the primary.
     * You can also set the primary flag explicitly:
     */

    int* primary = rebx_add_param(&sim->particles[0], "primary", REBX_TYPE_INT);
    *primary = 1;

    double tmax = 1000.;
    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}
