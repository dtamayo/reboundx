/**
 * Saving and loading simulations
 * 
 * This example demonstrates how to restart a simulation with all REBOUNDx effects and parameters.
 */
#include <stdio.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "core.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation(); // make a simple sim with star and 1 planet
    struct reb_particle p = {0};
    p.m     = 1.;   
    reb_add(sim, p);
    double m = 0.;
    double a = 1.e-4;
    double e = 0.2;
    double inc = 0.;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle p1 = reb_tools_orbit_to_particle(sim->G, p, m, a, e, inc, Omega, omega, f);
    reb_add(sim, p1);
    sim->dt = 1.e-8;
    sim->integrator = REB_INTEGRATOR_WHFAST;
    
    struct rebx_extras* rebx = rebx_attach(sim);
   
    // Add some (arbitrary) parameters to effects and particles
    
    struct rebx_force* gr = rebx_load_force(rebx, "gr");
    rebx_set_param_double(rebx, &gr->ap, "c", 3e4);
    rebx_set_param_int(rebx, &gr->ap, "gr_source", 0);
    rebx_add_force(rebx, gr);
    
    struct rebx_operator* mm = rebx_load_operator(rebx, "modify_mass");
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tau_mass", -1.e4);
    rebx_add_operator(rebx, mm);
    
    struct reb_particle* pptr = &sim->particles[1];
    rebx_set_param_double(rebx, &pptr->ap, "tau_e", 1.5);   // no meaning, just to test storage
    rebx_set_param_int(rebx, &pptr->ap, "max_iterations", 42);
    
    double* c = rebx_get_param(rebx, gr->ap, "c");
    int* gr_source = rebx_get_param(rebx, gr->ap, "gr_source");
    double* tau_mass = rebx_get_param(rebx, sim->particles[0].ap, "tau_mass");
    double* tau_e = rebx_get_param(rebx, sim->particles[1].ap, "tau_e");
    int* max_iterations = rebx_get_param(rebx, sim->particles[1].ap, "max_iterations");
    
    struct rebx_operator* integforce = rebx_load_operator(rebx, "integrate_force");
    rebx_set_param_pointer(rebx, &integforce->ap, "force", gr);
    rebx_add_operator(rebx, integforce);
    
    printf("c: Original = %f\n", *c);
    printf("gr_source: Original = %d\n", *gr_source);
    printf("tau_mass: Original = %f\n", *tau_mass);
    printf("tau_e: Original = %f\n", *tau_e);
    printf("max_iterations: Original = %d\n", *max_iterations);
    
    // We now have to save both a REBOUND binary (for the simulation) and a REBOUNDx one (for parameters and effects)
    reb_integrate(sim, 1.e-4);
    reb_output_binary(sim, "reb.bin");
    rebx_output_binary(rebx, "rebx.bin");
    
    rebx_free(rebx);
    reb_free_simulation(sim);
   
    // We now reload the simulation and the rebx instance (which adds previously loaded effects to the simulation)
    sim = reb_create_simulation_from_binary("reb.bin");
    rebx = rebx_create_extras_from_binary(sim, "rebx.bin");

    gr = rebx_get_force(rebx, "gr");
    mm = rebx_get_operator(rebx, "modify_mass");

    c = rebx_get_param(rebx, gr->ap, "c");
    gr_source = rebx_get_param(rebx, gr->ap, "gr_source");
    tau_mass = rebx_get_param(rebx, sim->particles[0].ap, "tau_mass");
    tau_e = rebx_get_param(rebx, sim->particles[1].ap, "tau_e");
    max_iterations = rebx_get_param(rebx, sim->particles[1].ap, "max_iterations");
     
    printf("c: Loaded = %f\n", *c);
    printf("gr_source: Loaded = %d\n", *gr_source);
    printf("tau_mass: Loaded = %f\n", *tau_mass);
    printf("tau_e: Loaded = %f\n", *tau_e);
    printf("max_iterations: Loaded = %d\n", *max_iterations);
    
    // You would now integrate as usual
    double tmax = 1.e-4;
    reb_integrate(sim, tmax);
    rebx_free(rebx);
    reb_free_simulation(sim);
}
