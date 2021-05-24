/**
 * Adding Post-Newtonian correction from general relativity as an operator.
 * 
 * This example shows how to add post-newtonian corrections to REBOUND simulations as an operator
 * (see REBOUNDx paper and corresponding IntegrateForce.ipynb jupyter notebook for more details).
 * If you have GLUT installed for the visualization, press 'w' and/or 'c' for a clearer view of
 * the whole orbit.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

double E0;
void heartbeat(struct reb_simulation* sim);

int main(int argc, char* argv[]){
    // We first set up a system similar to the EPIC system discussed in the REBOUNDx paper
    struct reb_simulation* sim = reb_create_simulation();
    sim->G = 4*M_PI*M_PI;

    struct reb_particle star = {0};
    star.m     = 0.93;   
    reb_add(sim, star);

    double m = 1.35e-6;
    double a = 0.013; // put planet close to enhance precession so it's visible in visualization (this would put planet inside the Sun!)
    double e = 0.01;
    double inc = 0.;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle planet1 = reb_tools_orbit_to_particle(sim->G, star, m, a, e, inc, Omega, omega, f);
    
    m = 1.2e-5;
    a = 0.107; // put planet close to enhance precession so it's visible in visualization (this would put planet inside the Sun!)
    e = 0.01;
    
    struct reb_particle planet2= reb_tools_orbit_to_particle(sim->G, star, m, a, e, inc, Omega, omega, f);
    reb_add(sim, planet1);
    reb_add(sim, planet2);
    reb_move_to_com(sim);
    
    sim->dt = 1.e-4;
    sim->heartbeat = heartbeat;
    sim->integrator = REB_INTEGRATOR_WHFAST;

    // Now we add GR. This is a velocity dependent force. With WHFast this would cause errors on long timescales, so we integrate the force in a separate step
    // See the REBOUNDx paper and the corresponding example in ipython_examples for more details
    // By adding a separate force step for GR and integrating across it we keep an oscillatory energy error
    
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* gr = rebx_load_force(rebx, "gr");
    rebx_set_param_double(rebx, &gr->ap, "c", 63197.8); // AU/yr

    struct rebx_operator* integrateforce = rebx_load_operator(rebx, "integrate_force");
    rebx_set_param_pointer(rebx, &integrateforce->ap, "force", gr);
    rebx_set_param_int(rebx, &integrateforce->ap, "integrator", REBX_INTEGRATOR_RK2);
    rebx_add_operator(rebx, integrateforce);

    double tmax = 10.;
    E0 = rebx_gr_hamiltonian(rebx, gr);
    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 0.1)){
        struct rebx_force* gr = rebx_get_force(sim->extras, "gr");
        double E = rebx_gr_hamiltonian(sim->extras, gr);
        printf("t=%f\tEnergy Error=%e\n", sim->t, fabs((E-E0)/E0));
    }
}
