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


void test_fragmenting_collisions(int type, int min_frag_test){
    // This function tests mass and momentum conservation for various setups.
    // **IMPORTANT** Tests for catching different collision regimes work with min_frag_mass = 0.05;
    struct reb_simulation* sim = reb_simulation_create(); //creates simulation
    sim->integrator = REB_INTEGRATOR_MERCURIUS;
    sim->collision = REB_COLLISION_DIRECT;
    sim->dt = 1;
    sim->rand_seed = 1;
    struct reb_particle star = {0};
    int eb_flag = 0; //eb_flag is to meant to flag where we expect an elastic bounce. This is useful to check IDs.

    // Add particles 
    // In Cases 
    switch (type){
        // In Cases 0 to 2, setup is to cause a merging event.
        case 0:
            reb_simulation_add_fmt(sim, "m r", 1.1, 1.0); // primary (slightly heavier)
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 1.0, 1.0, 2.5, -1.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            break;
        case 1: // order swapped
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 1.0, 1.0, 2.5, -1.0, 0.001, 0.001); 
            reb_simulation_add_fmt(sim, "m r", 1.1, 1.0); // primary (slightly heavier)
            break;
        case 2: // equal mass
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 1.0, 1.0, 2.5, -1.0, 0.001, 0.001); 
            reb_simulation_add_fmt(sim, "m r", 1.0, 1.0);
            break;
        // In cases 3 to 5, setup is to cause an erosion event.
        case 3: // Erosion, Non grazing
            reb_simulation_add_fmt(sim, "m r", 0.25, 1.0); // primary (slightly heavier)
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 0.20, 1.0, 10.0, -30.0, 0.001, 0.001); 
            break;
        case 4: // Erosion, Non grazing, order swapped
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 0.20, 1.0, 10.0, -30.0, 0.001, 0.001); 
            reb_simulation_add_fmt(sim, "m r", 0.25, 1.0); // primary (slightly heavier)
            break;
        case 5: // Erosion, Non grazing, equal mass
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 0.20, 1.0, 10.0, -30.0, 0.001, 0.001); 
            reb_simulation_add_fmt(sim, "m r", 0.20, 1.0);
            break;
        // In cases 6 to 11, setup is to cause a grazing event
        case 6:
            reb_simulation_add_fmt(sim, "m r", 0.25, 100.0); // primary (slightly heavier)
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.20, 100.0, 200.0, 150.0, -100.0, 0.001, 0.001); 
            break;
        case 7: // order swapped
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.20, 100.0, 200.0, 150.0, -100.0, 0.001, 0.001); 
            reb_simulation_add_fmt(sim, "m r", 0.25, 100.0); // primary (slightly heavier)
            break;
        case 8: // equal mass
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.20, 100.0, 200.0, 150.0, -100.0, 0.001, 0.001);
            reb_simulation_add_fmt(sim, "m r", 0.20, 100.0);
            break;
        case 9:
            reb_simulation_add_fmt(sim, "m r", 0.25, 100.0); //primary heavier, higher velocity
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.20, 100.0, 200.0, 150.0, -800.0, 0.001, 0.001);
            break;
        case 10:
            reb_simulation_add_fmt(sim, "m r vx", 0.25, 100.0, 100.0); // primary heavier, angle closer to 90deg
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.20, 100.0, 200.0, 50.0, -200.0, 0.001, 0.001);
            break;
        case 11:
            reb_simulation_add_fmt(sim, "m r", 0.25, 100.0); // primary heavier, angle closer to 90deg, smaller relative velocity 
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.20, 100.0, 200.0, 50.0, -100.0, 0.001, 0.001); 
            break;
        // In Cases 12 to 25, we're trying to catch different collision outcomes.
        // Refer to the flowchart (add reference) for specifications of different cases
        // Note: I have refered to cases with numbers, but to avoid confusion I will probably change to letters
        case 12: //merging (Case 1 (A))
            reb_simulation_add_fmt(sim, "m r", 0.501, 0.5); 
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.5, 0.5, 1.0, 0.9, -1.0, 0.001, 0.001); 
            break;
        case 13: //Elastic bounce (Case 8 (H))
            eb_flag = 1; //Since we expect an elastic bounce, eb_flag is 1.
            reb_simulation_add_fmt(sim, "m r", 0.101, 0.5); //mass ratio close to 1
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.1, 0.5, 1.0, 0.9, -1.0, 0.001, 0.001); 
            break;
        case 14: //Elastic bounce (Case 8 (H))
            eb_flag = 1;
            reb_simulation_add_fmt(sim, "m r", 0.4, 0.5); // 1:4 mass ratio, low velocity
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.1, 0.5, 1.0, 0.9, -1.0, 0.001, 0.001); 
            break;
        case 15: //Elastic bounce (Case 8 (H))
            eb_flag = 1;
            reb_simulation_add_fmt(sim, "m r", 0.4, 0.5); // 4:1 mass ratio, high velocity
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.1, 0.5, 8.0, 0.8, -10.0, 0.001, 0.001); 
            break;
        case 16: //Elastic bounce (Case 8 (H))
            eb_flag = 1;
            reb_simulation_add_fmt(sim, "m r", 100.0, 50.0); // 1:4 mass ratio, bigger masses 
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 25.0, 50.0, 75.0, 60.0, -80.0, 0.001, 0.001); 
            break;
        case 17: //Elastic bounce (Case 8 (H))
            eb_flag = 1;
            reb_simulation_add_fmt(sim, "m r", 0.1, 0.5); // 1:4 mass ratio, smaller mass
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.025, 0.5, 1.0, 0.7, -1.0, 0.001, 0.001); 
            break;
        case 18: //Elastic bounce (Case 8 (H))
            eb_flag = 1;
            reb_simulation_add_fmt(sim, "m r", 10.0, 0.5); // 10:1 mass ratio
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 1.0, 0.5, 8.0, 0.8, -10.0, 0.001, 0.001); 
            break;
        case 19: //Elastic bounce (Case 7 (G))
            eb_flag = 1;
            reb_simulation_add_fmt(sim, "m r", 0.1, 0.5); // 1:4 mass ratio, smaller mass, low velocity
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.025, 0.5, 1.0, 0.9, -1.0, 0.001, 0.001); 
            break;
        case 20: //Elastic bounce (Case 7 (G))
            eb_flag = 1;
            reb_simulation_add_fmt(sim, "m r", 0.4, 0.5); // 1:4 mass ratio, higher velocity
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.1, 0.5, 8.0, 0.9, -10.0, 0.001, 0.001); 
            break;
        case 21: //Elastic bounce (Case 5 (E))
            eb_flag = 1;
            reb_simulation_add_fmt(sim, "m r", 0.1, 0.5); // 1:10 mass ratio
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.01, 0.5, 50.0, 0.51, -60.0, 0.001, 0.001); 
            break;
        case 22: //Merging [M_rem too small] (Case 3 (C))
            reb_simulation_add_fmt(sim, "m r", 0.4, 0.5); // 1:4 mass ratio, 
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.1, 0.5, 8.0, 0.0, -10.0, 0.001, 0.001); 
            break;
        case 23: //Merging (Case 2 (B))
            reb_simulation_add_fmt(sim, "m r", 0.101, 0.5); // v_imp < v_crit
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.1, 0.5, 1.0, 0.9, -0.60, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            break;
        case 24: //Non grazing erosion (Case 4 (D))
            reb_simulation_add_fmt(sim, "m r", 4.0, 0.5); // 1:4 mass ratio, 
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 1.0, 0.5, 50.0, 0.0, -80.0, 0.001, 0.001); 
            break;
        case 25: //Grazing Erosion (Case 6 (F))
            reb_simulation_add_fmt(sim, "m r", 0.1, 0.5); // 1:10 mass ratio, higher velocity
            reb_simulation_add_fmt(sim, "m r x y vx vy vz", 0.01, 0.5, 200.0, 0.51, -250.0, 0.001, 0.001);
            break;
        //Cases 26 to 29 are edge cases
        case 26:
            //collision with the sun
            star.m = 1.00;
            star.r = 0.1; 
            reb_simulation_add(sim, star);
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 1e-7, 1e-5, 1.0, -2.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            break;
        case 27:
            //One object with zero mass
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 0.0, 1.0, 1.0, -1.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            reb_simulation_add_fmt(sim, "m r", 1.0, 1.0); // primary (slightly heavier)
            break;
        case 28:
            //The other object with zero mass
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 1.0, 1.0, 1.0, -1.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            reb_simulation_add_fmt(sim, "m r", 0.0, 1.0); // primary (slightly heavier)
            break;
        case 29:
            //Both objects with zero mass
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 0.0, 1.0, 1.0, -1.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            reb_simulation_add_fmt(sim, "m r", 0.0, 1.0);
            break;
    }

    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_collision_resolve* fragmenting = rebx_load_collision_resolve(rebx, "fragmenting_collisions");
    rebx_add_collision_resolve(rebx, fragmenting);
    //Test what happens if user doesn't define a min_frag_mass
    //Should be assigning default value
    if (min_frag_test == 0){
    }
    //Test what happens if user chooses a min frag mass
    else if(min_frag_test == 1){
        rebx_set_param_double(rebx, &fragmenting->ap, "fc_min_frag_mass", 0.01);
    }
    //Test when user sets a 0 min frag mass
    //Should be using default value and ignoring zero
    else if(min_frag_test == 2){
        rebx_set_param_double(rebx, &fragmenting->ap, "fc_min_frag_mass", 0.0);
    }

    // Assign all particles an initial id
    for(int i=0; i<sim->N; i++){
        rebx_fragmenting_collisions_set_new_id(sim, fragmenting, &sim->particles[i]);
    }
    
    struct reb_particle com_i = reb_simulation_com(sim); //initial center of mass
    reb_simulation_integrate(sim, 1);
    struct reb_particle com_f = reb_simulation_com(sim); //final center of mass

    //So far, we have 25 tests to check for conservation of mass and momentum
    //Since we don't want to check this for edge cases, I have excluded them here.
    int n_last_physical = 25;

    if(type <= n_last_physical){
        //printf("com_i - com_f = %e \n", (com_i.m-com_f.m)/com_i.m);
        assert(fabs((com_i.m-com_f.m)/com_i.m)<1e-10); // Check mass conservation 
        assert(fabs((com_i.vx-com_f.vx)/com_i.vx)<1e-10); // Check x momentum conservation 
        assert(fabs((com_i.vy-com_f.vy)/com_i.vy)<1e-10); // Check y momentum conservation 
        assert(fabs((com_i.vz-com_f.vz)/com_i.vz)<1e-10); // Check z momentum conservation 

        // ID checks
        int* fc_id_max = (int*) rebx_get_param(rebx, fragmenting->ap, "fc_id_max");
        assert(fc_id_max); // Make sure max ID has been assigned.
        assert(*fc_id_max != 0); // Make sure max ID is not 0.
        for(int i=0; i<sim->N; i++){
            int* fc_id = (int*) rebx_get_param(rebx, sim->particles[i].ap, "fc_id");
            assert(fc_id); // Make sure ID has been assigned.
            if(eb_flag == 0){ //If eb_flag = 1, then we have an elastic bounce and one of IDs will be zero, so no need for this test.
                assert(*fc_id != 0); // Make sure ID is not 0.
            }
            for(int j=0; j<sim->N; j++){
                if (i!=j){
                    int* fc_id2 = (int*) rebx_get_param(rebx, sim->particles[j].ap, "fc_id");
                    assert(fc_id2); // Make sure ID has been assigned.
                    assert(*fc_id != *fc_id2); // Make sure ID is unique
                }
            }
        }
    }
    reb_simulation_free(sim);
}


int main(int argc, char* argv[]) {
    for (int type=0;type<25;type++){
        for(int min_frag_test_flag=0;min_frag_test_flag<3;min_frag_test_flag++){
            test_fragmenting_collisions(type, min_frag_test_flag);
            printf("test_fragmenting_collisions(%d), round (%d) passed.\n", type, min_frag_test_flag);
        }
    }
}
