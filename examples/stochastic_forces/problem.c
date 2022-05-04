/**
 * Stochastic forces on a single planet
 * 
 * This example shows how to add a stochastic force to a single planet. 
 * As a result, the planet will undergo a random walk in its orbital
 * parameters. See also [Rein (2010)](https://ui.adsabs.harvard.edu/abs/2010PhDT.......196R/abstract) and [Rein and Papaloizou (2009)](https://ui.adsabs.harvard.edu/abs/2009A%26A...497..595R/abstract)
 *
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
    // Let's first setup the simulation.
    // Note that we are using the WHFast integrator with a fixed timestep. 
    // It's important to point out that the IAS15 integrator is not well 
    // suited for stochastic forces because it automatically reduces the 
    // timestep if it doesn't achieve an accuracy near machine precision. 
    // Because the stochastic forces are random, it might never converge.

    struct reb_simulation* sim = reb_create_simulation();
    sim->integrator     = REB_INTEGRATOR_WHFAST;
    sim->dt             = 1e-2;          // At ~100 AU, orbital periods are ~1000 yrs, so here we use ~1% of that, in sec
    sim->heartbeat      = heartbeat;
    
    reb_add_fmt(sim, "m", 1.); // Sun 
    reb_add_fmt(sim, "m a", 1e-3, 1.0); // Jupiter mass planet at 1AU
    reb_move_to_com(sim);

    // Next, we add the `stochastic_forces` module in REBOUNDx
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* sto = rebx_load_force(rebx, "stochastic_forces");
    rebx_add_force(rebx, sto);

    // We can now turn on stochastic forces for a particle by setting 
    // the field `kappa` to a finite value. This parameter determines 
    // the strength of the stochastic forces and is relative to the 
    // gravitational force that the particle experiences from the 
    // center of mass interior to its orbit. 
    rebx_set_param_double(rebx, &sim->particles[1].ap, "kappa", 1e-5);

    // The auto-correlation time of the stochastic forces is by default
    // the orbital period. We can set it to a fraction or multiple of
    // the orbital period by chanhging the `tau_kappa` parameter.

    // rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_kappa", 1.); // Not needed. Already the default.


    // Let's integrate the system.
    reb_integrate(sim, tmax);

    rebx_free(rebx);                                
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Periodically output the semi-major axis of the planet.
    if(reb_output_check(sim, 1.e2*M_PI*2.0)){
        reb_output_timing(sim, tmax);
        FILE* f = fopen("orbit.txt", "a+");
        struct reb_orbit o =  reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
        fprintf(f, "%.8e\t%.8e\n", sim->t, o.a);
        fclose(f);
    }
}
