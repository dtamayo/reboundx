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
void simple_drag(struct reb_simulation* const sim, struct rebx_effect* const effect)
{
    double* c = rebx_get_param(effect, "c");    // get parameters we want user to set
    double dt = sim->dt_last_done;              // important to check dt for integrators like IAS15 with adaptive timesteps

    if(c != NULL){                              // check to make sure c is set to avoid segmentation fault if not!
        sim->particles[2].vx *= 1. - (*c)*dt;
        sim->particles[2].vy *= 1. - (*c)*dt;
        sim->particles[2].vz *= 1. - (*c)*dt;
    }
}

// for a force we have to update the particle accelerations
void stark_force(struct reb_simulation* const sim, struct rebx_effect* const effect){
    double* c = rebx_get_param(effect, "c");    // get parameters we want user to set

    if(c != NULL){                              
        sim->particles[2].ax += (*c);           // make sure you += not =, which would overwrite other accelerations
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

    /* We now add our custom post_timestep_modification
     * We pass rebx, a name, and the function that should be called.
     * For a custom force, we also have to pass force_is_velocity_dependent,
     * which should be 1 if our function uses particle velocities, and 0 otherwise. 
     */

    struct rebx_effect* effect = rebx_add_custom_post_timestep_modification(rebx, "simple_drag", simple_drag);
    //struct rebx_effect* effect = rebx_add_custom_force(rebx, "stark_force", stark_force, 0);
    
    double* c = rebx_add_param(effect, "c", REBX_TYPE_DOUBLE);
    *c = 1.e-5;                                 // we wrote both our functions to read a parameter c, so we set it.
    
    double tmax = 5.e4;
    reb_integrate(sim, tmax);
    rebx_free(rebx);                            // Free all the memory allocated by rebx
}

