/**
 * Post-Newtonian correction from general relativity
 * 
 * This example shows how to add post-newtonian corrections to REBOUND simulations with reboundx.
 * If you have GLUT installed for the visualization, press 'w' and/or 'c' for a clearer view of
 * the whole orbit.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* r);
double E0;
struct rebx_effect* gr_params;
double tmax;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    struct rebx_extras* rebx = rebx_init(sim);
    sim->heartbeat        = heartbeat;
    struct reb_particle star = {0};
    star.m     = 1.;   
    star.hash  = reb_hash("star");
    reb_add(sim, star);

    double m = 1.e-3;
    double a = 1.;
    double e = 0.1;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle planet = reb_tools_orbit2d_to_particle(sim->G, star, m, a, e, omega, f);
    planet.hash = reb_hash("planet");
    reb_add(sim, planet);
    
    a = 10.;
    planet = reb_tools_orbit2d_to_particle(sim->G, star, m, a, e, omega, f);
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    // Could also add "gr" or "gr_full" here.  See documentation for details.
    gr_params = rebx_add(rebx, "gr");
   
    // Have to set speed of light in right units (set by G & initial conditions).  Here we use default units of AU/(yr/2pi)
    double* c = rebx_add_param(gr_params, "c", REBX_TYPE_DOUBLE);  
    *c = 1.e3;

    //int* noit = rebx_add_param(gr_params, "no_iteration", REBX_TYPE_INT);
    
    tmax = 10.;
    E0 = rebx_gr_hamiltonian(sim, gr_params);
    //reb_step(sim);
    reb_integrate(sim, tmax);
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    if (reb_output_check(sim, tmax/10.)){
        double E = rebx_gr_hamiltonian(sim, gr_params);
        struct reb_particle com = reb_get_com(sim);
        printf("%e\t%e\t%e\t%e\n",sim->t, fabs((E-E0)/E0), com.x, com.y);
        //printf("%e\t%e\t%e\n", E, reb_tools_energy(sim), fabs((E-reb_tools_energy(sim))/E));
    }
}
