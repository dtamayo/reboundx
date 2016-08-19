/**
 * Tracking a particle's minimum distance from the central star.
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    struct rebx_extras* rebx = rebx_init(sim);

    struct reb_particle star = {0};
    star.m     = 1.;   
    star.hash  = reb_hash("star");
    reb_add(sim, star);

    double m = 0.;
    double a = 1.;
    double e = 0.5;
    double omega = 0.;
    double f = M_PI;
    
    struct reb_particle planet = reb_tools_orbit2d_to_particle(sim->G, star, m, a, e, omega, f);
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    struct rebx_effect* track_min_distance = rebx_add(rebx, "track_min_distance");
  
    // Wee add a min_distance parameter to the particle whose distance we want to track, and set it
    // to a particular value. In any timestep that the distance drops below this value, the value is updated. 
    double* min_distance = rebx_add_param(&sim->particles[1], "min_distance", REBX_TYPE_DOUBLE);
    *min_distance = 5.;

    // By default distance is measured from sim->particles[0].  We can specify a different particle by a hash (unnecessary here):
    
    uint32_t* min_distance_from = rebx_add_param(&sim->particles[1], "min_distance_from", REBX_TYPE_UINT32);
    *min_distance_from = sim->particles[0].hash;

    double tmax = 10.;
    reb_integrate(sim, tmax);

    // At any point in the integration, we can check the `min_distance` parameter and output it as needed.
    min_distance = rebx_get_param(&sim->particles[1], "min_distance");
    printf("Particle's minimum distance from the star over the integration = %e\n", *min_distance);

    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}
