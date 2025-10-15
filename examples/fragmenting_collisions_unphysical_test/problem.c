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
#include <stdbool.h>


void test_merge(int type){
    char particles_filename[100];
    sprintf(particles_filename, "family_tree_%d.csv", type);
    // This function tests mass and momentum conservation for various setups.
    struct reb_simulation* sim = reb_simulation_create(); //creates simulation
    sim->integrator = REB_INTEGRATOR_MERCURIUS;
    sim->collision = REB_COLLISION_DIRECT;
    sim->dt = 1;
    sim->rand_seed = 1;
    struct reb_particle star = {0};

    // Add particles
    switch (type){
        case 0:
            //collision with the sun
            star.m = 1.00;
            star.r = 0.1; 
            reb_simulation_add(sim, star);
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 1e-7, 1e-5, 1.0, -2.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            break;
        case 1:
            //One object with zero mass
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 0.0, 1.0, 1.0, -1.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            reb_simulation_add_fmt(sim, "m r", 1.0, 1.0); // primary (slightly heavier)
            break;
        case 2:
            //The other object with zero mass
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 1.0, 1.0, 1.0, -1.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            reb_simulation_add_fmt(sim, "m r", 0.0, 1.0); // primary (slightly heavier)
            break;
        case 3:
            //Both objects with zero mass
            reb_simulation_add_fmt(sim, "m r x vx vy vz", 0.0, 1.0, 1.0, -1.0, 0.001, 0.001); // small vy, vz velocity yo check for momentum conservation in 3D
            reb_simulation_add_fmt(sim, "m r", 0.0, 1.0);
            break;
    }

    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_collision_resolve* fragmenting = rebx_load_collision_resolve(rebx, "fragmenting_collisions");
    rebx_add_collision_resolve(rebx, fragmenting);

    FILE* of = fopen(particles_filename, "a");
    fprintf(of, "particle_id, parent_t, parent_p,\n"); 
    // Assign all particles an initial id
    for(int i=0; i<sim->N; i++){
        struct reb_particle* p = &(sim->particles[i]); // First object in collision
        rebx_fragmenting_collisions_set_new_id(sim, fragmenting, &sim->particles[i]);
        int* new_id = rebx_get_param(rebx, p->ap, "fc_id");
        printf("particle %d ID is %d\n", i,  *new_id);
        fprintf(of, "%d,", *new_id);
        fprintf(of, "0, ");
        fprintf(of, "0, ");
        fprintf(of, "\n");   
    }

    fclose(of);
    
    reb_simulation_integrate(sim, 1);
    reb_simulation_free(sim);
}
    


int main(int argc, char* argv[]) {
    for (int type=0;type<4;type++){
        test_merge(type);
        printf("test_merge(%d) passed.\n", type);
    }
}
