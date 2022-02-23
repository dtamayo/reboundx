/**
 * Inner disk edge.
 * 
 * This example shows how to add an inner disk edge when using modify_orbits_forces or modify_orbits_direct.
 * See detailed annotations and explanations for the parameters in the InnerDiskEdge ipython example and 
 * Implemented Effects documentation for REBOUNDx.
 */
#include <stdio.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "core.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation(); 
    /*sim units are ('yr', 'AU', 'Msun')*/
    sim->G = 4*M_PI*M_PI;

    struct reb_particle star = {0};
    star.m     = 1.;  
    reb_add(sim, star);

    double m = 0.00001;
    double a = 1;
    double e = 0;
    double inc = 0.;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;

    struct reb_particle p1 = reb_tools_orbit_to_particle(sim->G, star, m, a, e, inc, Omega, omega, f);
    reb_add(sim, p1);
    reb_move_to_com(sim);

    sim->dt = 0.002;  //The period at the inner disk edge divided by 20, for a disk edge location at 0.1 AU
    sim->integrator = REB_INTEGRATOR_WHFAST;

    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* mof = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_set_param_double(rebx, &mof->ap, "ide_position", 0.1); 
    rebx_set_param_double(rebx, &mof->ap, "ide_width", 0.02); // Planet will stop within 0.02 AU of the inner disk edge
    rebx_add_force(rebx, mof);

    double tmax = 1.e4; 
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -tmax);        
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_e", -tmax/100.);     

    reb_integrate(sim, tmax);
    rebx_free(rebx);
    reb_free_simulation(sim);
}
