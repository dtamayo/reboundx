/**
 * Inner disk edge.
 * 
 * This example shows how to add an inner disk edge.
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

    sim->dt = 0.002;  //The period at inner disk edge divided by 20, for a disk edge location at 0.1 AU
    sim->integrator = REB_INTEGRATOR_WHFAST;

    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* inner_edge = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_set_param_double(rebx, &inner_edge->ap, "inner_disk_edge_position", 0.1); 
    rebx_set_param_double(rebx, &inner_edge->ap, "disk_edge_width", 0.02); //Calculated using a scale height value of 0.03. See Pichierri et al. 2018 for the equation
    rebx_add_force(rebx, inner_edge);

    double tmax = 1.e4; 
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -tmax);        
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_e", -tmax/100.);     

    reb_integrate(sim, tmax);
    rebx_free(rebx);
    reb_free_simulation(sim);
}