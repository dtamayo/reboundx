/**
 * Drag interaction from tides.
 *
 * This example shows how to add drag interaction due to slowly rotating tides raised on the primary body.
 * See also the corresponding Python example for more details.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* sim);
double tmax = 100.*2.*M_PI;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    sim->heartbeat = heartbeat;
    sim->integrator = REB_INTEGRATOR_IAS15;

    struct reb_particle star = {0};
    star.m     = 1.;
    reb_add(sim, star);

    struct reb_particle planet = {0};  // add a planet on a circular orbit (with default units where G=1)
    planet.m = 1e-6;
    planet.x = 1.;
    planet.vy = 1.;
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    struct rebx_extras* rebx = rebx_attach(sim);  // first initialize rebx
    struct rebx_force* tides = rebx_load_force(rebx, "tides_drag"); // add our new force
    // rebx_set_param_double(rebx, &sim->particles[1].ap, "migration_tau", 1000); // set particle parameter
    rebx_add_force(rebx, tides);

    // Overwrite planet output file w/ col headers
    system("rm -f planet.txt"); // remove existing file
    FILE* file = fopen("planet.txt","a");
    fprintf(file, "Time(yrs)\t\tMass(Msun)\t\tSemi-major Axis(AU)\t\tEccentricity\t\tInclination(Radians)\t\tLongitude_of_Ascending_Node(Radians)\t\tArgument_of_Periapsis(Radians))\t\tTrue_Anomaly(Radians)\n");
    fclose(file);
    
    reb_integrate(sim, tmax);
}

void heartbeat(struct reb_simulation* sim){
    if (reb_output_check(sim, tmax/1000)){
        // retrieve Sun particle
        struct reb_particle sun = sim->particles[0];
        // retrieve Earth particle
        struct reb_particle ep = sim->particles[1];
        struct reb_orbit eo  = reb_tools_particle_to_orbit(sim->G, ep, sun);
        double t = sim->t;
        double m = ep.m;
        double a = eo.a;
        double e = eo.e;
        double inc = eo.inc;
        double Omega = eo.Omega;
        double omega = eo.omega;
        double f = eo.f;
        FILE* file = fopen("planet.txt","a");

        reb_output_timing(sim, tmax);
        reb_integrator_synchronize(sim);
        fprintf(file,"%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n",t,m,a,e,inc,Omega,omega,f);
        fclose(file);
    }
}
