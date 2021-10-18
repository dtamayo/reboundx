/**
 * Radiation forces on circumplanetary dust
 * 
 * This example shows how to integrate circumplanetary
 * dust particles under the action of radiation forces using IAS15.
 * We use Saturn's Phoebe ring as an example, a distant ring of debris, 
 * The output is custom, outputting the semi-major axis of 
 * every dust particle relative to the planet. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* r);

double tmax = 1e10;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Setup constants 
    sim->integrator     = REB_INTEGRATOR_IAS15;
    sim->G              = 6.674e-11;    // Use SI units
    sim->dt             = 1e4;          // Initial timestep in sec
    sim->N_active       = 2;            // Only the sun and the planet affect other particles gravitationally
    sim->heartbeat      = heartbeat;
    sim->usleep     = 1000;             // Slow down integration (for visualization only)
    
    // sun
    struct reb_particle sun = {0};
    sun.m  = 1.99e30;                   // mass of Sun in kg
    reb_add(sim, sun);
    
    // Saturn (simulation set up in Saturn's orbital plane, i.e., inc=0, so only need 4 orbital elements 
    double mass_sat = 5.68e26;          // Mass of Saturn 
    double a_sat = 1.43e12;             // Semimajor axis of Saturn in m
    double e_sat = 0.056;               // Eccentricity of Saturn
    double pomega_sat = 0.;             // Angle from x axis to pericenter
    double f_sat = 0.;                  // True anomaly of Saturn
    double inc = 0.;
    double Omega = 0.;
    struct reb_particle saturn = reb_tools_orbit_to_particle(sim->G, sun, mass_sat, a_sat, e_sat, inc, Omega, pomega_sat, f_sat);

    reb_add(sim, saturn);

    // Add REBOUNDx
    struct rebx_extras* rebx = rebx_attach(sim); 
    struct rebx_force* rad = rebx_load_force(rebx, "radiation_forces");
    double c = 3.e8;                    // speed of light in SI units 
    rebx_set_param_double(rebx, &rad->ap, "c", c);
    
    // Will assume particles[0] is the radiation source by default. You can also add a flag to a particle explicitly
    rebx_set_param_int(rebx, &sim->particles[0].ap, "radiation_source", 1); 

    /* Dust particles
     Here we imagine particles launched from Saturn's irregular Satellite Phoebe.
     Such grains will inherit the moon's orbital elements (e.g. Tamayo et al. 2011) 
        
     In order for a particle to feel radiation forces, we have to set their beta parameter, 
     the ratio of the radiation pressure force to the gravitional force from the star (Burns et al. 1979). 
     We do this in two ways below.*/
        
    double a_dust = 1.30e10;            // semimajor axis of satellite Phoebe, in m
    double e_dust = 0.16;               // eccentricity of Phoebe
    double inc_dust = 175.*M_PI/180.;   // inclination of Phoebe to Saturn's orbital plane
    double Omega_dust = 0.;             // longitude of ascending node
    double omega_dust = 0.;             // argument of pericenter
    double f_dust = 0.;                 // true anomaly

    // We first set up the orbit and add the particles
    double m_dust = 0.;                 // treat dust particles as massless
    struct reb_particle p = reb_tools_orbit_to_particle(sim->G, sim->particles[1], m_dust, a_dust, e_dust, inc_dust, Omega_dust, omega_dust, f_dust); 
    reb_add(sim, p); 

    // For the first particle we simply specify beta directly.
    rebx_set_param_double(rebx, &sim->particles[2].ap, "beta", 0.1); 

    // We now add a 2nd particle on the same orbit, but set its beta using physical parameters.  
    struct reb_particle p2 = reb_tools_orbit_to_particle(sim->G, sim->particles[1], 0., a_dust, e_dust, inc_dust, Omega_dust, omega_dust, f_dust); 
    reb_add(sim, p2);

    /* REBOUNDx has a convenience function to calculate beta given the gravitational constant G, the star's luminosity and mass, and the grain's physical radius, density and radiation pressure coefficient Q_pr (Burns et al. 1979). */
   
    // Particle parameters
    double radius = 1.e-5;              // in meters
    double density = 1.e3;              // kg/m3 = 1g/cc 
    double Q_pr = 1.;                   // Equals 1 in limit where particle radius >> wavelength of radiation
    double L = 3.85e26;                 // Luminosity of the sun in Watts

    double beta = rebx_rad_calc_beta(sim->G, c, sim->particles[0].m, L, radius, density, Q_pr);
    rebx_set_param_double(rebx, &sim->particles[3].ap, "beta", beta); 

    printf("Particle 2 has beta = %f\n", beta);

    reb_move_to_com(sim);

    printf("Time\t\tEcc (p)\t\tEcc (p2)\n");
    reb_integrate(sim, tmax);
    rebx_free(rebx);                /* free memory allocated by REBOUNDx */
}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 1.e8)){
        const struct reb_particle* particles = sim->particles;
        const struct reb_particle saturn = particles[1];
        const double t = sim->t;
        struct reb_orbit orbit = reb_tools_particle_to_orbit(sim->G, sim->particles[2], saturn); /* calculate orbit of particles[2] around Saturn */
        double e2 = orbit.e;
        orbit = reb_tools_particle_to_orbit(sim->G, sim->particles[3], saturn); 
        double e3 = orbit.e;
        printf("%e\t%f\t%f\n", t, e2, e3);
    }
}
