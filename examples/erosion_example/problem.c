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


void test_erosion(int type){
    // This function tests mass and momentum conservation for various setups.
    struct reb_simulation* sim = reb_simulation_create(); //creates simulation
    sim->integrator = REB_INTEGRATOR_MERCURIUS;
    sim->collision = REB_COLLISION_DIRECT;
    sim->dt = 1;
    sim->rand_seed = 1;

    // Add particles
    switch (type){
        case 0:
            reb_simulation_add_fmt(sim, "m r", 0.25, 1.0); // primary (slightly heavier)
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 0.20, 1.0, 10.0, -30.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            break;
        case 1: // order swapped
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 0.20, 1.0, 10.0, -30.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            reb_simulation_add_fmt(sim, "m r", 0.25, 1.0); // primary (slightly heavier)
            break;
        case 2: // equal mass
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 0.20, 1.0, 10.0, -30.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            reb_simulation_add_fmt(sim, "m r", 0.20, 1.0);
            break;
    }

    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_collision_resolve* fragmenting = rebx_load_collision_resolve(rebx, "fragmenting_collisions");
    rebx_add_collision_resolve(rebx, fragmenting);
    
    reb_simulation_integrate(sim, 1);
}
    


int main(int argc, char* argv[]) {
    test_erosion(0);
}
