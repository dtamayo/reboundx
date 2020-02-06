/**
 * Sun-Earth Output
 *
 * This example integrates a simple Sun-Earth system.
 * Heartbeat function outputs Earth orbital data for
 * later import and and graphical analysis in Python.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* r);
double tmax;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->dt = 1./20.;              // 1/20 Earth's period in yrs
    tmax = 8e4;                  // 4 Myr
    r->G = 4*M_PI*M_PI;          // in AU^3 / Msun / yr^2.
    r->ri_whfast.safe_mode = 0;  // Turn off safe mode. Need to call reb_integrator_synchronize() before outputs. 
    r->ri_whfast.corrector = 11; // 11th order symplectic corrector
    r->integrator = REB_INTEGRATOR_WHFAST;
    r->heartbeat = heartbeat;
    r->exact_finish_time = 1;    // Finish exactly at tmax in reb_integrate(). Default is already 1.

    // Add Sun to sim
    struct reb_particle sun = {0}; // initialize w/ zeroes
    sun.m = 1;                     // in Msun
    reb_add(r, sun);
    // Add Earth to sim
    struct reb_orbit eo = {0};
    double e_mass = .000002988;
    eo.a = 1.0;                    // in AU
    struct reb_particle ep = reb_tools_orbit_to_particle(r->G, sun, e_mass, eo.a, eo.e, eo.inc, eo.Omega, eo.omega, eo.f);
    reb_add(r, ep);

    // Overwrite planet output file w/ col headers
    system("rm -f planet.txt"); // remove existing file
    FILE* file = fopen("planet.txt","a");
    fprintf(file, "Time(yrs)\t\tMass(Msun)\t\tSemi-major Axis(AU)\t\tEccentricity\t\tInclination(Radians)\t\tLongitude_of_Ascending_Node(Radians)\t\tArgument_of_Periapsis(Radians))\t\tTrue_Anomaly(Radians)\n");
    fclose(file);

    // Overwrite COM output file
    system("rm -f COM.txt"); // remove existing file
    FILE* com_file = fopen("COM.txt","a");
    fprintf(com_file, "Time(yrs)\t\tx(AU)\t\ty(AU)\t\tz(AU)\n");
    fclose(com_file);

    // Run simulation
    reb_move_to_com(r);
    reb_integrate(r, tmax);
}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 1000.)){
        // retrieve COM "particle"
        struct reb_particle com = reb_get_com(r);
        // retrieve Sun particle
        struct reb_particle sun = r->particles[0];
        // retrieve Earth particle
        struct reb_particle ep = r->particles[1];
        struct reb_orbit eo  = reb_tools_particle_to_orbit(r->G, ep, sun);
        double t = r->t;
        double m = ep.m;
        double a = eo.a;
        double e = eo.e;
        double inc = eo.inc;
        double Omega = eo.Omega;
        double omega = eo.omega;
        double f = eo.f;
        FILE* file = fopen("planet.txt","a");
        FILE* com_file = fopen("COM.txt","a");

        reb_output_timing(r, tmax);
        reb_integrator_synchronize(r);
        fprintf(file,"%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n",t,m,a,e,inc,Omega,omega,f);
        fclose(file);
        fprintf(com_file,"%e\t\t%e\t\t%e\t\t%e\t\t\n",t,com.x,com.y,com.z);
        fclose(com_file);
    }
}

