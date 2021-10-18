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
    struct rebx_extras* rebx = rebx_attach(sim);

    struct reb_particle star = {0};
    star.m     = 1.;   
    star.hash  = reb_hash("star");
    reb_add(sim, star);

    double m = 0.;
    double a = 1.;
    double e = 0.9;
    double inc = 0.;
    double Omega = 0.;
    double omega = 0.;
    double f = M_PI;
    
    struct reb_particle planet = reb_tools_orbit_to_particle(sim->G, star, m, a, e, inc, Omega, omega, f);
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    struct rebx_operator* track_min_distance = rebx_load_operator(rebx, "track_min_distance");
    rebx_add_operator(rebx, track_min_distance);
  
    // Wee add a min_distance parameter to the particle whose distance we want to track, and set it
    // to a particular value. In any timestep that the distance drops below this value, the value is updated. 
    rebx_set_param_double(rebx, &sim->particles[1].ap, "min_distance", 5.);

    // By default distance is measured from sim->particles[0].  We can specify a different particle by a hash (unnecessary here):
    
    rebx_set_param_uint32(rebx, &sim->particles[1].ap, "min_distance_from", sim->particles[0].hash);


    struct reb_orbit orbit = {0};
    rebx_set_param_pointer(rebx, &sim->particles[1].ap, "min_distance_orbit", &orbit);

    double tmax = 10.;
    reb_integrate(sim, tmax);

    // At any point in the integration, we can check the `min_distance` parameter and output it as needed.
    double* min_distance = rebx_get_param(rebx, sim->particles[1].ap, "min_distance");
    struct reb_orbit* orbitptr = rebx_get_param(rebx, sim->particles[1].ap, "min_distance_orbit");

    printf("Particle's minimum distance from the star over the integration = %e\n", *min_distance);
    printf("Semimajor axis and eccentricity at closest approach: a=%.2f, e=%.2f\n", orbitptr->a, orbitptr->e);

    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}
