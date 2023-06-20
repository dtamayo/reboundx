/**
 * Adding Lense-Thirring effect
 * 
 * This example shows how to add the Lense-Thirring effect to a simulation.
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
    sim->G = 4*M_PI*M_PI; // units of AU, yr, Msun
    struct reb_particle star = {0};
    star.m     = 1.;
    reb_add(sim, star);
    double omega = 90.361036076; //solar rotation rate in rad/year
    double C_I = 0.06884; //solar moment of inertia prefactor
    double R_eq = 0.00465247264; //solar equatorial radius in AU

    struct reb_particle planet = {0};  // add a planet on a circular orbit (with default units where G=1)
    const double mp = 1.7e-7;   // approximate values for mercury in units of Msun and AU
    const double a = 0.39;
    const double e = 0.21;
    reb_add_fmt(sim, "m a e", mp, a, e);

    struct rebx_extras* rebx = rebx_attach(sim);  // first initialize rebx
    struct rebx_force* lense  = rebx_load_force(rebx, "lense_thirring"); // add our new force
    rebx_add_force(rebx, lense);
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.x=0, .y=0, .z=omega});
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", C_I*star.m*R_eq*R_eq);

    // Have to set speed of light in right units (set by G & initial conditions).  Here we use default units of AU/(yr/2pi)
    rebx_set_param_double(rebx, &lense->ap, "lt_c", 63241.077); // speed of light in AU/yr

    double tmax = 1000.;
    reb_integrate(sim, tmax);
    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    printf("Pericenter precession rate = %.3f arcsec/century\n", orb.pomega*206265/10.);
}
