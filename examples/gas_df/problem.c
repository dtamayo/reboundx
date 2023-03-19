/**
 * Gas dynamical friction
 * 
 * This example shows how to add the gas_df force.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    sim->integrator = REB_INTEGRATOR_BS;

    struct reb_particle bh = {0};
    bh.m     = 4e6;   
    reb_add(sim, bh);

    double m = 1;
    double a = 206000; 
    double e = 0.01;
    double inc = 0.17;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle star = reb_tools_orbit_to_particle(sim->G, bh, m, a, e, inc, Omega, omega, f);
    reb_add(sim, star);
    reb_move_to_com(sim);
    
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* gdf = rebx_load_force(rebx, "gas_df");
    rebx_add_force(rebx, gdf);

    rebx_set_param_double(rebx, &gdf->ap, "gas_df_rhog", 0.20);    
    rebx_set_param_double(rebx, &gdf->ap, "gas_df_alpha_rhog", -1.5);
    rebx_set_param_double(rebx, &gdf->ap, "gas_df_cs", 20);
    rebx_set_param_double(rebx, &gdf->ap, "gas_df_alpha_cs", -0.5);
    rebx_set_param_double(rebx, &gdf->ap, "gas_df_xmin", 0.045);
    rebx_set_param_double(rebx, &gdf->ap, "gas_df_hr", 0.01);
    rebx_set_param_double(rebx, &gdf->ap, "gas_df_Qd", 5.0);

    double delta_t = 6.28e5;
    for (int i = 0; i < 100; i++){
        reb_integrate(sim, sim->t + delta_t);
        struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
        printf("%f %f %f %e\n", sim->t, o.a, o.e, o.inc);
    }

    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}
