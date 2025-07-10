/**
 * Merging collisions
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
    sim->integrator = REB_INTEGRATOR_LEAPFROG;
    sim->collision = REB_COLLISION_DIRECT;
    sim->dt = 1.0;

    // Add particle at origin:
    struct reb_particle p = {0};
    p.r = 1.;
    reb_simulation_add(sim, p);

    // Shift particle and add a copy:
    p.r = 1.;
    p.x = 2.5;
    p.vx = -1;
    reb_simulation_add(sim, p);

    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_collision_resolve* fragmenting = rebx_load_collision_resolve(rebx, "fragmenting_collisions");
    rebx_add_collision_resolve(rebx, fragmenting);

    reb_simulation_integrate(sim, 1); //integrates system for 1 timestep
    
    printf("Number of particles left in simulation: %d \n", sim->N); // Only one particle left.
    
}
