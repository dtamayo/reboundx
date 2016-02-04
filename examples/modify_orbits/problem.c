/**
 * Migration and other orbit modifications
 *
 * This example shows how to add migration, eccentricity damping
 * and pericenter precession to a REBOUND simulation.  If you have
 * GLUT installed for visualization, press 'w' to see the orbits
 * as wires.  You can zoom out by holding shift, holding down the mouse
 * and dragging.  Press 'c' to better see migration/e-damping.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Setup constants
    sim->dt             = 0.012;        // initial timestep.
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

    // There are two options for how to modify orbits.  You would only choose one (comment the other out).  
    // modify_orbits_forces doesn't have precession implemented yet.

    rebx_add_modify_orbits_direct(rebx);    // directly update particles' orbital elements each timestep
    //rebx_add_modify_orbits_forces(rebx);  // add forces that orbit-average to give exponential a and e damping

    // Set the timescales for each particle.  Parameter getter and setter functions always take the address of the particle (&)
    
    rebx_set_tau_a(&sim->particles[1], -1.e5); // add semimajor axis damping on inner planet (e-folding timescale)
    rebx_set_tau_omega(&sim->particles[1], -1.e4); // add linear precession on inner planet (precession period) 
    rebx_set_tau_e(&sim->particles[2], -1.e4); // add eccentricity damping on particles[2] (e-folding timescale)

    printf("Semimajor axis damping timescale for inner planet is %f.\n", -1.*rebx_get_tau_a(&sim->particles[1]));
    printf("Precession timescale for inner planet is %f.\n", -1.*rebx_get_tau_omega(&sim->particles[1]));
    printf("Eccentricity damping timescale for outer planet is %f.\n", -1.*rebx_get_tau_e(&sim->particles[2]));
   
    /* One can also adjust a coupling parameter between eccentricity and semimajor axis damping.  We use the parameter p
     * as defined by Deck & Batygin (2015).  The default p=0 corresponds to no coupling, while p=1 corresponds to e-damping
     * at constant angular momentum.  This is currently only implemented for modify_orbits_direct (not modify_orbits_forces).
     *
     * Additionally, the damping by default is done in Jacobi coordinates.  If you'd prefer to use barycentric or heliocentric
     * coordinates, set rebx->modify_orbits_forces.coordinates to BARYCENTRIC or HELIOCENTRIC, respectively. (also works for
     * modify_orbits_direct)*/

    struct rebx_params_modify_orbits* params = rebx_get_effect_params(rebx, MODIFY_ORBITS_DIRECT);
    params->p = 1.;
    params->coordinates = HELIOCENTRIC;

    double tmax = 5.e4;
    reb_integrate(sim, tmax);
    rebx_free(rebx);    // Free all the memory allocated by rebx
}
