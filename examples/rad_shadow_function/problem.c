#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]) {
    struct reb_simulation* sim = reb_simulation_create();
    
    sim->gravity = REB_GRAVITY_NONE; 
    sim->G = 6.674e-11;  
    sim->dt = 1.0;       

    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* rf = rebx_load_force(rebx, "radiation_forces");
    
    double c_light = 299792458.0;
    rebx_set_param_double(rebx, &rf->ap, "c", c_light);
    
    int shadow_model = 2; // 2 = Setting up the conical model
    rebx_set_param_int(rebx, &rf->ap, "shadow_model", shadow_model);

    rebx_add_force(rebx, rf);
    // Radiation Source in the Centre 
    struct reb_particle sun = {0};
    sun.m = 1.989e30;
    sun.x = 0.0; sun.y = 0.0; sun.z = 0.0;
    sun.vx = 0.0; sun.vy = 0.0; sun.vz = 0.0;
    reb_simulation_add(sim, sun);
    
    rebx_set_param_int(rebx, &sim->particles[0].ap, "radiation_source", 1);
    double r_sun = 6.957e8;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "R_eq", r_sun);

    // Occulting body (Earth) at 1 AU
    double au = 1.496e11;
    struct reb_particle earth = {0};
    earth.m = 5.972e24;
    earth.x = au; earth.y = 0.0; earth.z = 0.0;
    earth.vx = 0.0; earth.vy = 0.0; earth.vz = 0.0;
    reb_simulation_add(sim, earth);
    
    int shadow_creator = 1; // Setting the particle up as a occulting body
    rebx_set_param_int(rebx, &sim->particles[1].ap, "shadow_creator", shadow_creator);
    double r_earth = 6.371e6;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "R_eq", r_earth);

    struct reb_particle pc = {0};
    pc.m = 0.0;
    pc.x = au + 10000e3; pc.y = 0.0; pc.z = 0.0;
    pc.vx = 0.0; pc.vy = 0.0; pc.vz = 0.0;
    reb_simulation_add(sim, pc);
    
    double beta = 0.1;
    rebx_set_param_double(rebx, &sim->particles[2].ap, "beta", beta);

    struct reb_particle pd = {0};
    pd.m = 0.0;
    pd.x = au - 10000e3; pd.y = 0.0; pd.z = 0.0;
    pd.vx = 0.0; pd.vy = 0.0; pd.vz = 0.0;
    reb_simulation_add(sim, pd);
    
    rebx_set_param_double(rebx, &sim->particles[3].ap, "beta", beta);


    reb_simulation_step(sim);

    // Verification 1: Total umbra
    printf("Acceleration in UMBRA (Expected: 0.000000e+00) : %e\n", sim->particles[2].ax);

    // Verification 2: No shadow
    double dist_light = sim->particles[3].x;
    double expected_a_light = beta * sim->G * sim->particles[0].m / (dist_light * dist_light);
    printf("Acceleration in LIGHT (Expected: %e) : %e\n\n", expected_a_light, sim->particles[3].ax);

    reb_simulation_free(sim);
    reb_simulation_free(rebx);
    return 0;
}