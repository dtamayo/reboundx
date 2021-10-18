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

/* There are two qualitatively different options for adding in effects. One can add a force, which evaluates accelerations to add
 * on top of gravitational accelerations every timestep and will be integrated numerically.
 * Alternatively one can apply operators before and/or after each REBOUND timestep that update particle properties
 * (masses, positions/velocities, or other parameters).
 * 
 * Here's an example of a force, for whichi we have to update the particle accelerations in the passed particles array (not sim->particles!)
 * Function must have the same prototype as below.
 * We can have the user assign parameters to the force structure in main(), and our force function read and/or update those parameters each timestep when called.
 */

void stark_force(struct reb_simulation* const sim, struct rebx_force* const starkforce, struct reb_particle* const particles, const int N){
    double* starkconst = rebx_get_param(sim->extras, starkforce->ap, "starkconst");    // get parameters we want user to set

    if(starkconst != NULL){                              
        particles[1].ax += (*starkconst);           // make sure you += not =, which would overwrite other accelerations
    }
}

/* Here's a very simple example of an operator to change the planet's orbit
 *
 * For a post_timestep_modification we update the particle states (positions, velocities, masses etc.)
 * Function must have the same prototype as below, passing simulation and operator pointers (which can be used to hold parameters), 
 * and the timestep over which the operator should act.
 */

void simple_drag(struct reb_simulation* const sim, struct rebx_operator* const dragoperator, const double dt){
    double* dragconst = rebx_get_param(sim->extras, dragoperator->ap, "dragconst");    // get parameters we want user to set

    if(dragconst != NULL){                              
        sim->particles[1].vx *= 1. - (*dragconst)*dt;
        sim->particles[1].vy *= 1. - (*dragconst)*dt;
        sim->particles[1].vz *= 1. - (*dragconst)*dt;
    }
}


int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    struct reb_particle p = {0}; 
    p.m     = 1.;   
    reb_add(sim, p); 

    double m = 0.;
    double a1 = 1.;
    double e = 0.4;
    double inc = 0.;
    double Omega = 0.;
    double omega = 0.;
    double f = 0.;

    struct reb_particle p1 = reb_tools_orbit_to_particle(sim->G, p, m, a1, e, inc, Omega, omega, f);
    reb_add(sim,p1);
    reb_move_to_com(sim);

    struct rebx_extras* rebx = rebx_attach(sim);  // first initialize rebx

    /* We now add our custom functions. We do this in the same way as we add REBOUNDx forces and operators.
     * We'll still get back a force and operator struct, respectively, with a warning that the name wasn't found 
     * in REBOUNDx and that we're reponsible for setting their type flags and function pointers.
     */

    struct rebx_force* stark = rebx_create_force(rebx, "stark_force");

    /* We first set a flag for whether our force only depends on particle positions (REBX_FORCE_POS) 
     * or whether it depends on velocities (or velocities and positions, REBX_FORCE_VEL)
     */
    stark->force_type = REBX_FORCE_VEL;
    stark->update_accelerations = stark_force;  // set the function pointer to what we wrote above
    rebx_add_force(rebx, stark);                // Now it's initialized, add to REBOUNDx
    
    struct rebx_operator* drag = rebx_create_operator(rebx, "simple_drag"); // Now we create our custom operator
   
    /* We first set the operator_type enum for whether our operator modifies dynamical variables (positions, velocities
     * or masses, REBX_OPERATOR_UPDATER), or whether it's just passively recording the state of the simulation, or 
     * updating parameters that don't feed back on the dynamics (REBX_OPERATOR_RECORDER)
     */

    drag->operator_type = REBX_OPERATOR_UPDATER;
    drag->step_function = simple_drag;  // set function pointer to what we wrote above
    rebx_add_operator(rebx, drag);      // Now it's initialized, add to REBOUNDx

    /* Now we set the parameters that our custom functions above need
     * Before setting them, we need to register them with their corresponding types
     */
   
    rebx_register_param(rebx, "starkconst", REBX_TYPE_DOUBLE);
    rebx_register_param(rebx, "dragconst", REBX_TYPE_DOUBLE);

    rebx_set_param_double(rebx, &stark->ap, "starkconst", 1.e-5);
    rebx_set_param_double(rebx, &drag->ap, "dragconst", 1.e-5);
   
    /* To simplify things for other users, we can always add the new effects and register parameters in the REBOUNDx
     * repo itself. Feel free to send a pull request or contact me (tamayo.daniel@gmail.com) about adding it
     */

    double tmax = 5.e4;
    reb_integrate(sim, tmax);
    rebx_free(rebx);                            // Free all the memory allocated by rebx
}
