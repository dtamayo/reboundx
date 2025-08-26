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


void test_merge(int type){
    // This function tests mass and momentum conservation for various setups.
    struct reb_simulation* sim = reb_simulation_create(); //creates simulation
    sim->integrator = REB_INTEGRATOR_MERCURIUS;
    sim->collision = REB_COLLISION_DIRECT;
    sim->dt = 1;
    sim->rand_seed = 1;

    // Add particles
    switch (type){
        case 0:
            reb_simulation_add_fmt(sim, "m r", 1.1, 1.0); // primary (slightly heavier)
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 1.0, 1.0, 2.5, -1.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            break;
        case 1: // order swapped
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 1.0, 1.0, 2.5, -1.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            reb_simulation_add_fmt(sim, "m r", 1.1, 1.0); // primary (slightly heavier)
        case 2: // equal mass
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 1.0, 1.0, 2.5, -1.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            reb_simulation_add_fmt(sim, "m r", 1.0, 1.0); // primary (slightly heavier)
            break;
    }

    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_collision_resolve* fragmenting = rebx_load_collision_resolve(rebx, "fragmenting_collisions");
    rebx_add_collision_resolve(rebx, fragmenting);
    
    struct reb_particle com_i = reb_simulation_com(sim); //initial center of mass
    reb_simulation_integrate(sim, 1);
    struct reb_particle com_f = reb_simulation_com(sim); //final center of mass

    assert(sim->N == 1); // Check that merge did occur
    assert(fabs((com_i.m-com_f.m)/com_i.m)<1e-16); // Check mass conservation 
    assert(fabs((com_i.vx-com_f.vx)/com_i.vx)<1e-16); // Check x momentum conservation 
    assert(fabs((com_i.vy-com_f.vy)/com_i.vy)<1e-16); // Check y momentum conservation 
    assert(fabs((com_i.vz-com_f.vz)/com_i.vz)<1e-16); // Check z momentum conservation 
}
    


int main(int argc, char* argv[]) {
    for (int type=0;type<3;type++){
        test_merge(0);
        printf("test_merge(%d) passed.\n", type);
    }
}
