/**
 * Highly eccentric orbits
 *
 * This example uses the IAS15 integrator to simulate
 * a very eccentric planetary orbit. The integrator
 * automatically adjusts the timestep so that the pericentre passages
 * resolved with high accuracy.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

double tmax;
void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->G            = 1;        // Gravitational constant
    r->integrator        = REB_INTEGRATOR_IAS15;
    r->heartbeat        = heartbeat;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->gravity = REB_GRAVITY_NONE;
    r->usleep = 10.;

    struct reb_particle tp = {0};
    tp.x = 1.5;
    tp.vy = 0.4;//1/sqrt(1.5); // for circular orbit

    reb_add(r, tp); 

    struct rebx_extras* rebx = rebx_attach(r);
    // Could also add "gr" or "gr_full" here.  See documentation for details.
    struct rebx_force* ephem_forces = rebx_load_force(rebx, "ephemeris_forces");
    rebx_add_force(rebx, ephem_forces); 
    // Have to set speed of light in right units (set by G & initial conditions).  Here we use default units of AU/(yr/2pi)
    rebx_set_param_int(rebx, &ephem_forces->ap, "N_ephem", 2);

    tmax            = 10000.;

    reb_integrate(r, tmax);
    r->dt = 1.e-4;
    reb_step(r);
    reb_step(r);
    reb_step(r);

}

void heartbeat(struct reb_simulation* r){
    if(reb_output_check(r,tmax/10000.)){        // outputs to the screen
        reb_output_timing(r, tmax);
    }
}

