/*
 * Adding J2 to Phoebe.
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_simulation_create();

    struct reb_particle star = {0};
    star.m     = 1.;   
    star.hash  = reb_hash("star");
    reb_simulation_add(sim, star);

    double m = 0.;
    double a = 1.; 
    double e = 0.2;
    double inc = 0.;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle planet = reb_particle_from_orbit(sim->G, star, m, a, e, inc, Omega, omega, f);
    planet.hash = reb_hash("planet");
    reb_simulation_add(sim, planet);
    reb_simulation_move_to_com(sim);
    
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* gh = rebx_load_force(rebx, "j2");
    rebx_add_force(rebx, gh);

    rebx_set_param_double(rebx, &sim->particles[0].ap, "J2", 0.1);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "R_eq", 0.01);
   
    double tmax = 1.e4;
    reb_simulation_integrate(sim, tmax); 
    rebx_free(rebx);
    reb_simulation_free(sim);
}

