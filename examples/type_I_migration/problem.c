/**
 * Type I migration.
 * 
 * This example shows how to add Type I migration.
 */
#include <stdio.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "core.h"

//DO I need the heartbeat void function?

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation(); 
    sim.units = ('yr', 'AU', 'Msun') //how do I write this in c?

    struct reb_particle star = {0};  //should I explicitly add a star?
    star.m     = 1.;   
    reb_add(sim, star);

    struct reb_particle p = {0};
    p.m     = 1.;  // what is one here? Is this the star included??
    reb_add(sim, p);

    double m = 0.00001;
    double a1 = 0.5;
    double a2 = 0.85;
    double e = 0;
    double inc = 0.;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;

    struct reb_particle p1 = reb_tools_orbit_to_particle(sim->G, p, m, a1, e, inc, Omega, omega, f);
    struct reb_particle p2 = reb_tools_orbit_to_particle(sim->G, p, m, a2, e, inc, Omega, omega, f);
    reb_add(sim, p1);
    reb_add(sim, p2);
    reb_move_to_com(sim);

    sim->dt = 0.002;  //The period at inner disk edge divided by 20, for a disk edge location at 0.1 AU
    sim->integrator = REB_INTEGRATOR_WHFAST;

    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* type_I_mig = rebx_load_force(rebx, "type_I_migration");
    rebx_set_param_double(rebx, &force->ap, "inner_disk_edge_position", 0.1); //DO I need to say that it is a double?
    rebx_set_param_double(rebx, &force->ap, "disk_edge_width", 0.02); //calculated using the scale height value given below
    rebx_set_param_double(rebx, &force->ap, "flaring_index", 0.25);
    rebx_set_param_double(rebx, &force->ap, "surface_density_exponent", 1);
    rebx_set_param_double(rebx, &force->ap, "initial_surface_density", 0.00011255); //in units of Msun/AU^2 which is roughly 1000g/cm^2
    rebx_set_param_double(rebx, &force->ap, "scale_height", 0.03);
    rebx_add_force(rebx, type_I_mig);

    double tmax = 1.e4;  // note that one can calculate the timescale of semi-major axis of the outer planet then set the integration time to twice this value 

    reb_integrate(sim, tmax);
    rebx_free(rebx);
    reb_free_simulation(sim);
}