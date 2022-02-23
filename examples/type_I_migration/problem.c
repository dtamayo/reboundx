/**
 * Type I migration.
 * 
 * This example shows how to add Type I migration.
 * See detailed annotations and explanations for the parameters in the TypeIMigration ipython example and 
 * Implemented Effects documentation for REBOUNDx.
 */
#include <stdio.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "core.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation(); 
    /* sim units are ('yr', 'AU', 'Msun') */
    sim->G = 4*M_PI*M_PI;

    struct reb_particle star = {0};  
    star.m     = 1.;   
    reb_add(sim, star);

    double m = 0.00001;
    double a1 = 0.5;
    double a2 = 0.85;
    double e = 0;
    double inc = 0.;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;

    struct reb_particle p1 = reb_tools_orbit_to_particle(sim->G, star, m, a1, e, inc, Omega, omega, f);
    struct reb_particle p2 = reb_tools_orbit_to_particle(sim->G, star, m, a2, e, inc, Omega, omega, f);
    reb_add(sim, p1);
    reb_add(sim, p2);
    reb_move_to_com(sim);

    sim->dt = 0.002;  //The period at inner disk edge divided by 20, for a disk edge location at 0.1 AU
    sim->integrator = REB_INTEGRATOR_WHFAST;

    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* type_I_mig = rebx_load_force(rebx, "type_I_migration");
    rebx_set_param_double(rebx, &type_I_mig->ap, "ide_position", 0.3); 
    rebx_set_param_double(rebx, &type_I_mig->ap, "ide_width", 0.02); //Calculated using the scale height value given below
    rebx_set_param_double(rebx, &type_I_mig->ap, "tIm_flaring_index", 0.25);
    rebx_set_param_double(rebx, &type_I_mig->ap, "tIm_surface_density_exponent", 1);
    rebx_set_param_double(rebx, &type_I_mig->ap, "tIm_surface_density_1", 0.00011255); //In units of solarmass/(AU^2) which is roughly 1000g/cm^2
    rebx_set_param_double(rebx, &type_I_mig->ap, "tIm_scale_height_1", 0.03);
    rebx_add_force(rebx, type_I_mig);

    double tmax = 1.e4;  // Note that one can calculate the timescale of the semi-major axis of the outer planet then set the integration time to twice this value 

    reb_integrate(sim, tmax);
    rebx_free(rebx);
    reb_free_simulation(sim);
}
