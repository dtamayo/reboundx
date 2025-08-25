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
    sim->rand_seed = 1;
    double tot_mass_i = 0;
    
    // Add particle at origin:
    struct reb_particle p = {0};
    p.r = 1.;
    p.m = 1.;
    printf("Particle 0 mass = %f \n", p.m);
    printf("Particle 0 xvelocity = %f \n", p.vx);
    printf("Particle 0 x = %f \n", p.x);
    tot_mass_i += p.m;
    reb_simulation_add(sim, p);

    // Shift particle and add a copy:
    p.r = 1.;
    p.x = 2.5;
    p.vx = -1;
    p.m = 1.;
    printf("Particle 1 mass = %f \n", p.m);
    printf("Particle 1 xvelocity = %f \n", p.vx);
    printf("Particle 1 x = %f \n", p.x);
    tot_mass_i += p.m;
    reb_simulation_add(sim, p);

    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_collision_resolve* merging = rebx_load_collision_resolve(rebx, "merging_collisions");
    rebx_add_collision_resolve(rebx, merging);

    reb_simulation_integrate(sim, 1); //integrates system for 1 timestep
    
    printf("Number of particles left in simulation: %d \n", sim->N); // Only one particle left.
    struct reb_particle p_temp;
    double tot_mass_f = 0; //final total mass
    for (int i=0; i < sim->N ; i++){ //this should technically be i < N, but I'm printing more cause they exist?!
        p_temp = sim->particles[i];
        tot_mass_f += p_temp.m;
        printf("Particle %d mass = %f \n", i, p_temp.m);
        printf("Particle %d xvelocity = %f \n", i, p_temp.vx);
        printf("Particle %d x = %f \n", i, p_temp.x);
    }
    printf("Total mass after collision = %f \n", tot_mass_f);
    
}
