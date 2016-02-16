/**
 * Adding custom post-timestep modifications and forces.
 *
 * This allows the user to use the built-in functions of REBOUNDx
 * but also include their own specialised functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

// A very simple post timestep modification to change the planet's orbit
// The custom function *has* to take (struct reb_simulation* const sim, struct rebx_effect* const effect)
// (you can just ignore the effect struct in your function body if you don't want to use it)

void simple_drag(struct reb_simulation* const sim, struct rebx_effect* const effect);

// Each effect has its own parameters structure.  You can define your own:
struct custom_params{
    double coefficient;
};

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Setup constants
    sim->dt             = 0.012;    // initial timestep.
    sim->integrator = REB_INTEGRATOR_IAS15;

    struct reb_particle p = {0}; 
    p.m     = 1.;   
    reb_add(sim, p); 

    double m = 0.;
    double a1 = 1.;
    double a2 = 2.;
    double e = 0.4;
    double omega = 0.;
    double f = 0.;

    struct reb_particle p1 = reb_tools_orbit2d_to_particle(sim->G, p, m, a1, e, omega, f);
    struct reb_particle p2 = reb_tools_orbit2d_to_particle(sim->G, p, m, a2, e, omega, f);
    reb_add(sim,p1);
    reb_add(sim,p2);
    reb_move_to_com(sim);

    struct rebx_extras* rebx = rebx_init(sim);  // first initialize rebx

    struct custom_params* params = malloc(sizeof(*params));
    params->coefficient = 1.e-5;
    
    // We now add our custom post_timestep_modification
    // We pass rebx, the function that should be called, and a pointer to our params structure.
    // You don't have to use a params structure (you would just pass NULL for the last argument).
    // In this example, you could hardcode 1.e-5 in simple_drag and not use a params structure.
    rebx_add_custom_post_timestep_modification(rebx, simple_drag, params);
    
    double tmax = 5.e4;
    reb_integrate(sim, tmax);
    rebx_free(rebx);                // Free all the memory allocated by rebx
}

void simple_drag(struct reb_simulation* const sim, struct rebx_effect* const effect)
{
    struct custom_params* params = effect->paramsPtr;
    double c = params->coefficient;
 
    sim->particles[2].vx *= 1. - c*sim->dt_last_done;
    sim->particles[2].vy *= 1. - c*sim->dt_last_done;
    sim->particles[2].vz *= 1. - c*sim->dt_last_done;
}
