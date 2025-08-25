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
    sim->rand_seed = 1;
    double tot_mass_i = 0; //initial total mass

    // Add particle 0
    struct reb_particle p = {0};
    p.r = 1.;
    p.m = 0.2;
    tot_mass_i += p.m;
    printf("Particle 0 mass = %f \n", p.m);
    printf("Particle 0 xvelocity = %f \n", p.vx);
    printf("Particle 0 x = %f \n", p.x);
    reb_simulation_add(sim, p);


    // Add particle 1
    p.r = 1.;
    p.x = 10.;
    int deltax = p.x; //save this to get number of timesteps for integration
    //set p.vx = -30 for fragmentation to occur
    //set p.vx = -2 (or some other low number) for merging to occur
    p.vx = -10.;
    int v = p.vx; //save this to get n_timesteps for integration
    p.m = 0.2;
    tot_mass_i += p.m;
    printf("Particle 1 mass = %f \n", p.m);
    printf("Particle 1 xvelocity = %f \n", p.vx);
    printf("Particle 1 x = %f \n", p.x);
    reb_simulation_add(sim, p);


    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_collision_resolve* fragmenting = rebx_load_collision_resolve(rebx, "fragmenting_collisions");
    rebx_add_collision_resolve(rebx, fragmenting);

    printf("Number of particles before collision: %d \n", sim->N);
    printf("Total mass before collision = %f \n", tot_mass_i);

    //integrate system until collision happens
    int n_timesteps = (int)(fabs(deltax) / fabs(v) / (sim->dt));
    printf("n_tsteps = %d \n", n_timesteps);
    reb_simulation_integrate(sim, n_timesteps + 1);


    printf("Number of particles left in simulation: %d \n", sim->N);
    struct reb_particle p_temp;
    double tot_mass_f = 0; //final total mass
    for (int i=0; i <= sim->N + 2; i++){ //this should technically be i < N, but I'm printing more cause they exist?!
        p_temp = sim->particles[i];
        tot_mass_f += p_temp.m;
        printf("Particle %d mass = %f \n", i, p_temp.m);
        printf("Particle %d xvelocity = %f \n", i, p_temp.vx);
        printf("Particle %d x = %f \n", i, p_temp.x);
    }
    printf("Total mass after collision = %f \n", tot_mass_f);
    
}
