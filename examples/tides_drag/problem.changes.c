/**
 * Precession from tides.
 *
 * This example shows how to add precession due to tides raised on either the primary, the orbiting bodies, or both.
 * See also the corresponding Python example for more details.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* sim);
double tmax = 4e6;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    sim->heartbeat = heartbeat;
    sim->dt = 1./20.;              // 1/20 Earth's period in yrs
    sim->G = 4*M_PI*M_PI;          // in AU^3 / Msun / yr^2.

    // Add Sun to sim
    struct reb_particle sun = {0}; // initialize w/ zeroes
    sun.m = 1;                     // in Msun
    reb_add(sim, sun);

    // Add Earth to sim
    struct reb_orbit eo = {0};
    double e_mass = 2.988e-6;
    eo.a = 1.0;                    // in AU
    struct reb_particle ep = reb_tools_orbit_to_particle(sim->G, sun, e_mass, eo.a, eo.e, eo.inc, eo.Omega, eo.omega, eo.f);
    reb_add(sim, ep);

    reb_move_to_com(sim);
    
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* tides = rebx_load_force(rebx, "tides_precession");
    rebx_add_force(rebx, tides);
   
    // Have to set R_tides (physical radius) and k1 (apsidal motion constant, half the tidal Love number k2) parameters. 
    // Could just set one set to consider tides raised on one body only.

    // rebx_set_param_double(rebx, &sim->particles[0].ap, "R_tides", 0.005);
    // rebx_set_param_double(rebx, &sim->particles[0].ap, "k1", 0.03);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "R_tides", 7.e4/1.5e8); // radius in same units of distance (here G=1, so AU)
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k1", 0.3);

    /* By default, implementation assumes particles[0] is the primary.
     * You can also set the tides_primary flag explicitly (don't have to set it to a value):
     */

    // rebx_set_param_int(rebx, &sim->particles[0].ap, "tides_primary", 0);

    // Overwrite planet output file w/ col headers
    system("rm -f planet.txt"); // remove existing file
    FILE* file = fopen("planet.txt","a");
    fprintf(file, "Time(yrs)\t\tMass(Msun)\t\tSemi-major Axis(AU)\t\tEccentricity\t\tInclination(Radians)\t\tLongitude_of_Ascending_Node(Radians)\t\tArgument_of_Periapsis(Radians))\t\tTrue_Anomaly(Radians)\n");
    fclose(file);

    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
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