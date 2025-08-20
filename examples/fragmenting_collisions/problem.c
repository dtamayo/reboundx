/**
 * Fragmenting Collisions
 *
 * A simple example showing how to use a REBOUNDx collision module.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]) {

    struct reb_simulation* sim = reb_simulation_create(); //creates simulation
    sim->integrator = REB_INTEGRATOR_MERCURIUS;
    sim->collision = REB_COLLISION_DIRECT;
    sim->dt = 0.1;

    // Add particle at origin:
    struct reb_particle p = {0};
    p.r = 1.;
    p.m = 0.2;
    reb_simulation_add(sim, p);

    // Shift particle and add a copy:
    p.r = 1.;
    p.x = 10;
    p.vx = -30;
    p.m = 0.2;
    reb_simulation_add(sim, p);

    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_collision_resolve* fragmenting = rebx_load_collision_resolve(rebx, "fragmenting_collisions");
    rebx_add_collision_resolve(rebx, fragmenting);

    reb_simulation_integrate(sim, 30); //integrates system for 1 timestep
    
    printf("Number of particles left in simulation: %d \n", sim->N); // Only one particle left.
    struct reb_particle p1_f = sim->particles[0]; 
    struct reb_particle p2_f = sim->particles[1]; 
    struct reb_particle p3_f = sim->particles[2]; 
    if(sim->N == 3){
    printf("Mass of p1: %f \n", p1_f.m); // Only one particle left.
    printf("Mass of p2: %f \n", p2_f.m); // Only one particle left.
    printf("Mass of p3: %f \n", p3_f.m); // Only one particle left.
    }
    
}
