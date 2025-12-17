/**
 * Chaotic dynamical tides model
 *
 * Close pericenter passages excite fundamental (f) modes in the planet.
 * This effect models the chaotic growth in energy in these modes, and their non-linear dissipation at large enough amplitudes.
 * See the corresponding ipython example, as well as the paper, for more explanations along the way of the various parameters and assumptions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* r);
double tmax = 2.e3*2*M_PI; // years

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_simulation_create();
    sim->heartbeat = heartbeat;

    struct reb_particle sun = {0}; 
    sun.m = 1.;
    reb_simulation_add(sim, sun);

    double a0 = 1.5;               // in AU
    double e0 = 0.987;
    double Rp = 1.6*4.67e-4;       // 1.6 Jup radii in AU
    double Mp = 1.e-3;             // 1 Jup mass in Msun
    struct reb_particle p = reb_particle_from_orbit(sim->G, sun, Mp, a0, e0, 0, 0, 0, 0);
    reb_simulation_add(sim, p);
    sim->particles[1].r = Rp;      // We must set both the mass and radius of the planet
    reb_simulation_move_to_com(sim);
    
    // Add dynamical tides
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* tides = rebx_load_force(rebx, "tides_dynamical");
    rebx_add_force(rebx, tides);

    // See the ipython example for various options that can be set

    // Run simulation
    reb_simulation_integrate(sim, tmax);
    rebx_free(rebx); 
    reb_simulation_free(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Periodically output the semi-major axis and eccentricity of the planet.
    if(reb_simulation_output_check(sim, 1.e2*M_PI*2.0)){
        struct reb_orbit o =  reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
        printf("%.3f\t%.3f\t%.3f\n", sim->t, o.a, o.e);
    }
}

