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

/* A very simple post timestep modification to change the planet's orbit
 * All custom function *have* to take (struct reb_simulation* const sim, struct rebx_effect* const effect)
 * (you can just ignore the effect struct in your function body if you don't want to use it)
 */

// for a post_timestep_modification we update the particle states (positions, velocities, masses etc.)
// function must have the same prototype as below, passing simulation and effect pointers, the timestep to apply,
// and an enum telling you whether the function is being called before or after the timestep (not used here).
void simple_drag(struct reb_simulation* const sim, struct rebx_effect* const effect, const double dt, enum rebx_timing timing){
    double* c = rebx_get_param(effect, "c");    // get parameters we want user to set

    if(c != NULL){                              // check to make sure c is set to avoid segmentation fault if not!
        sim->particles[2].vx *= 1. - (*c)*dt;
        sim->particles[2].vy *= 1. - (*c)*dt;
        sim->particles[2].vz *= 1. - (*c)*dt;
    }
}

// for a force we have to update the particle accelerations in the passed particles array (not sim->particles!)
// function must have the same prototype as below
void stark_force(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N){
    double* b = rebx_get_param(effect, "b");    // get parameters we want user to set

    if(b != NULL){                              
        particles[2].ax += (*b);           // make sure you += not =, which would overwrite other accelerations
    }
}

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
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

    /* We now add our custom operator
     * We pass rebx, a name for the effect, and the function that should be called.
     * For a custom force, we also have to pass force_is_velocity_dependent,
     * which should be 1 if our function uses particle velocities, and 0 otherwise. 
     */

    struct rebx_effect* drag = rebx_add_custom_operator(rebx, "simple_drag", simple_drag);
    struct rebx_effect* stark = rebx_add_custom_force(rebx, "stark_force", stark_force, 0);
    
    double* c = rebx_add_param(drag, "c", REBX_TYPE_DOUBLE);
    double* b = rebx_add_param(stark, "b", REBX_TYPE_DOUBLE);
    *c = 1.e-5;                                 
    *b = 1.e-5;
    
    double tmax = 5.e4;
    reb_integrate(sim, tmax);
    rebx_free(rebx);                            // Free all the memory allocated by rebx
}
