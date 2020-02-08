/**
 * Drag interaction from tides.
 *
 * This example shows how to add drag interaction due 
 * to slowly rotating tides raised on the primary body.
 * In particular, this simulates post-main sequence
 * tidal interaction between the Earth and Sun near
 * its tip-RGB phase.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    double tmax = 7e4;    // in yr
    sim->dt = 1./20.;     // 1/20 Earth's period (year)
    sim->G = 4*M_PI*M_PI; // in AU^3 / Msun / yr^2.
    
    struct reb_particle sun = {0}; // initialize w/ zeroes
    sun.m = 0.86;                  // in Msun
    reb_add(sim, sun);

    struct reb_orbit eo = {0};     // for Earth
    eo.a = 1.0;                    // in AU
    double e_mass = 2.988e-6;
    struct reb_particle ep = reb_tools_orbit_to_particle(sim->G, sun, e_mass, eo.a, eo.e, eo.inc, eo.Omega, eo.omega, eo.f);
    reb_add(sim, ep);
    
    // Add REBOUNDx Additional Effect
    struct rebx_extras* rebx = rebx_attach(sim);                                // initialize rebx
    struct rebx_force* tides = rebx_load_force(rebx, "tides_drag");             // add force
    rebx_add_force(rebx, tides);

    // Set optional and required params
    rebx_set_param_int(rebx, &sim->particles[0].ap, "tides_primary", 1);        // opt
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tides_Omega", 0);       // opt
    rebx_set_param_double(rebx, &sim->particles[0].ap, "R_tides", 0.78);        // req'd
    rebx_set_param_double(rebx, &sim->particles[0].ap, "luminosity", 869.5);    // req'd
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tides_lambda2", 0.023); // req'd
       
    // Run simulation
    reb_move_to_com(sim);
    reb_integrate(sim, tmax);
}
