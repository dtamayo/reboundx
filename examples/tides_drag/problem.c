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
double tmax;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Setup constants
    sim->dt = 1./20.;              // 1/20 Earth's period (yrs)
    // sim->dt = 0.0120423;           // 1/20 Mercury's period (yrs)
    tmax = 8e4;                    // 80 Kyr
    sim->G = 4*M_PI*M_PI;          // in AU^3 / Msun / yr^2.
    sim->ri_whfast.safe_mode = 0;  // Turn off safe mode. Need to call reb_integrator_synchronize() before outputs. 
    sim->ri_whfast.corrector = 11; // 11th order symplectic corrector
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->heartbeat = heartbeat;

    // Add Sun to sim
    struct reb_particle sun = {0}; // initialize w/ zeroes
    sun.m = 0.86;                  // in Msun
    reb_add(sim, sun);
    // Add Earth to sim
    struct reb_orbit eo = {0};
    double e_mass = .000002988;
    eo.a = 1.0;                    // in AU
    struct reb_particle ep = reb_tools_orbit_to_particle(sim->G, sun, e_mass, eo.a, eo.e, eo.inc, eo.Omega, eo.omega, eo.f);
    reb_add(sim, ep);
    
    // Add REBOUNDx Additional Effect
    struct rebx_extras* rebx = rebx_attach(sim);  // first initialize rebx
    struct rebx_force* tides = rebx_load_force(rebx, "tides_drag"); // add our new force
    // struct rebx_force* tides = rebx_load_force(rebx, "migration_force"); // add our new force
    // rebx_set_param_double(rebx, &sim->particles[1].ap, "migration_tau", 8e4); // set particle parameter
    rebx_add_force(rebx, tides);

    // Overwrite planet output file w/ col headers
    system("rm -f planet.txt"); // remove existing file
    FILE* file = fopen("planet.txt","a");
    fprintf(file, "Time(yrs)\t\tMass(Msun)\t\tSemi-major Axis(AU)\t\tEccentricity\t\tInclination(Radians)\t\tLongitude_of_Ascending_Node(Radians)\t\tArgument_of_Periapsis(Radians))\t\tTrue_Anomaly(Radians)\n");
    fclose(file);
    
    reb_move_to_com(sim);
    reb_integrate(sim, tmax);
}

void heartbeat(struct reb_simulation* sim){
    if (reb_output_check(sim, 1000.)){
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
