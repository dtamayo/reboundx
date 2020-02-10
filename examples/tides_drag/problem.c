/**
 * Drag interaction from tides.
 *
 * This example shows how to add drag interaction due to slowly rotating tides raised on the primary body.
 * The effect assumes all affected particles follow co-planar, circular orbits with respect to the primary’s rotational pole.
 * In particular, this simulates post-main sequence tidal interactions between the Earth and Sun near its tip-RGB phase.
 * See also the corresponding Python example for more details.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* sim);
double tmax = 7e4; // in yr

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Define constants
    sim->dt = 1./20.;     // 1/20 Earth's period (yr)
    sim->G = 4*M_PI*M_PI; // in AU^3 / Msun / yr^2.
    sim->heartbeat = heartbeat;
    
    struct reb_particle sun = {0}; // initialize w/ zeroes
    sun.m = 0.86;                  // post-MS in Msun
    reb_add(sim, sun);

    struct reb_orbit eo = {0};     // for Earth
    eo.a = 1.0;                    // in AU
    double e_mass = 2.988e-6;      // if planets don't have mass, they won't be affected.
    struct reb_particle ep = reb_tools_orbit_to_particle(sim->G, sun, e_mass, eo.a, eo.e, eo.inc, eo.Omega, eo.omega, eo.f);
    reb_add(sim, ep);
    reb_move_to_com(sim);
    
    // Add REBOUNDx Additional Effect
    struct rebx_extras* rebx = rebx_attach(sim);                                // initialize rebx
    struct rebx_force* tides = rebx_load_force(rebx, "tides_drag");             // add force
    rebx_add_force(rebx, tides);

    // Have to set R_tides (physical radius), luminosity (in units consistent with G), and tides_lambda2 (a coefficient that depends on the properties of the primary’s convective envelope) parameters.
    rebx_set_param_double(rebx, &sim->particles[0].ap, "R_tides", 0.78);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "luminosity", 869.5);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tides_lambda2", 0.023);

    // By default, implementation assumes particles[0] is the primary.
    // You can also set the tides_primary flag explicitly (don't have to set it to a value):
    rebx_set_param_int(rebx, &sim->particles[0].ap, "tides_primary", 1);
    // By default, implementation assumes primary rotational angular velocity is 0 unless otherwise specified.
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tides_Omega", 0);

    // Run simulation
    reb_integrate(sim, tmax);
    rebx_free(rebx); // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    if (reb_output_check(sim, 1000.)){
        // Output progress every 1000 yrs
        reb_output_timing(sim, tmax);
    }
}
