/**
 * Stochastic forces on a single planetk
 * 
 * This example shows how to add a stochastic force to a single planet. 
 * As a result, the planet will undergo a random walk in its orbital
 * parameters.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* r);

double tmax = 1e4*M_PI*2.0; // 1e4 orbits

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Setup constants
    sim->integrator     = REB_INTEGRATOR_WHFAST;
    sim->dt             = 1e-2;          // At ~100 AU, orbital periods are ~1000 yrs, so here we use ~1% of that, in sec
    sim->heartbeat      = heartbeat;
    
    reb_add_fmt(sim, "m", 1.); // Sun 
    reb_add_fmt(sim, "m a", 1e-3, 1.0); // Jupiter mass planet at 1AU
    reb_move_to_com(sim);
    
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* sto = rebx_load_force(rebx, "stochastic_forces");
    rebx_add_force(rebx, sto);

    rebx_set_param_double(rebx, &sim->particles[1].ap, "D", 1e-5);

    reb_integrate(sim, tmax);

    rebx_free(rebx);                                
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 1.e2*M_PI*2.0)){
        reb_output_timing(sim, tmax);
        FILE* f = fopen("orbit.txt", "a+");
        struct reb_orbit o =  reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
        fprintf(f, "%.8e\t%.8e\n", sim->t, o.a);
        fclose(f);
    }
}
