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

    // Assign all particles an initial id
    for(int i=0; i<sim->N; i++){
        rebx_fragmenting_collisions_set_new_id(sim, fragmenting, &sim->particles[i]);
    }
    
    struct reb_particle com_i = reb_simulation_com(sim); //initial center of mass
    reb_simulation_integrate(sim, 1);
    struct reb_particle com_f = reb_simulation_com(sim); //final center of mass

    assert(fabs((com_i.m-com_f.m)/com_i.m)<1e-14); // Check mass conservation 
    assert(fabs((com_i.vx-com_f.vx)/com_i.vx)<1e-13); // Check x momentum conservation 
    assert(fabs((com_i.vy-com_f.vy)/com_i.vy)<1e-13); // Check y momentum conservation 
    assert(fabs((com_i.vz-com_f.vz)/com_i.vz)<1e-13); // Check z momentum conservation 
   
    // ID checks
    int* fc_id_max = (int*) rebx_get_param(rebx, fragmenting->ap, "fc_id_max");
    assert(fc_id_max); // Make sure max ID has been assigned.
    assert(*fc_id_max != 0); // Make sure max ID is not 0.
    for(int i=0; i<sim->N; i++){
        int* fc_id = (int*) rebx_get_param(rebx, sim->particles[i].ap, "fc_id");
        assert(fc_id); // Make sure ID has been assigned.
        assert(*fc_id != 0); // Make sure ID is not 0.
        for(int j=0; j<sim->N; j++){
            if (i!=j){
                int* fc_id2 = (int*) rebx_get_param(rebx, sim->particles[j].ap, "fc_id");
                assert(fc_id2); // Make sure ID has been assigned.
                assert(*fc_id != *fc_id2); // Make sure ID is unique
            }
        }
    }
    reb_simulation_free(sim);
}
    


int main(int argc, char* argv[]) {
    for (int type=0;type<3;type++){
        test_erosion(type);
        printf("test_erosion(%d) passed.\n", type);
    }
}
