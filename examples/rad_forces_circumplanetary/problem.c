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
    /* Setup constants */
    sim->integrator     = REB_INTEGRATOR_IAS15;
    sim->G              = 6.674e-11;    /* Use SI units */
    sim->dt             = 1e4;          /* Initial timestep in sec */
    sim->N_active       = 2;            /* Only the sun and the planet affect other particles gravitationally */
    sim->heartbeat      = heartbeat;
    sim->usleep     = 1000;             /* Slow down integration (for visualization only) */
    
    /* sun */
    struct reb_particle sun = {0};
    sun.m  = 1.99e30;                   /* mass of Sun in kg */
    reb_add(sim, sun);
    
    struct rebx_extras* rebx = rebx_init(sim); 
    double c = 3.e8;                    /* speed of light in SI units */
    int source_index = 0;               /* Index of the particle that is the source of radiation. */
    struct rebx_params_radiation_forces* params = rebx_add_radiation_forces(rebx, source_index, c); 
    
    /* Saturn (simulation set up in Saturn's orbital plane, i.e., inc=0, so only need 4 orbital elements */
    double mass_sat = 5.68e26;          /* Mass of Saturn */
    double a_sat = 1.43e12;             /* Semimajor axis of Saturn in m */
    double e_sat = 0.056;               /* Eccentricity of Saturn */
    double pomega_sat = 0.;             /* Angle from x axis to pericenter */
    double f_sat = 0.;                  /* True anomaly of Saturn */
    struct reb_particle saturn = reb_tools_orbit2d_to_particle(sim->G, sun, mass_sat, a_sat, e_sat, pomega_sat, f_sat);

    reb_add(sim, saturn);

    /* Dust particles
     Here we imagine particles launched from Saturn's irregular Satellite Phoebe.
     Such grains will inherit the moon's orbital elements (e.g. Tamayo et al. 2011) 
        
     In order for a particle to feel radiation forces, we have to set their beta parameter, 
     the ratio of the radiation pressure force to the gravitional force from the star (Burns et al. 1979). 
     We do this in two ways below.*/
        
    double a_dust = 1.30e10;/* semimajor axis of satellite Phoebe, in m */
    double e_dust = 0.16;   /* eccentricity of Phoebe */
    double inc_dust = 175.*M_PI/180.;   /* inclination of Phoebe to Saturn's orbital plane */
    double Omega_dust = 0.; /* longitude of ascending node */
    double omega_dust = 0.; /* argument of pericenter */
    double f_dust = 0.;     /* true anomaly */

    /* We first set up the orbit and add the particles*/
    double m_dust = 0.;     /* treat dust particles as massless */
    struct reb_particle p = reb_tools_orbit_to_particle(sim->G, sim->particles[1], m_dust, a_dust, e_dust, inc_dust, Omega_dust, omega_dust, f_dust); 
    reb_add(sim, p); 

    /* We can only call REBOUNDx parameter setters on particles in the sim->particles array, so we do this after adding the particle to sim. For the first particle we simply specify beta directly.*/
    double beta = 0.1;
    rebx_set_param_double(&sim->particles[2], "beta", beta);    

    /* We now add a 2nd particle on the same orbit, but set its beta using physical parameters.  Note that callingreb_add(sim, p) a second time will add a second particle on the same orbit, but their beta paramters will be linked (see read_me_first example) */ 
    struct reb_particle p2 = reb_tools_orbit_to_particle(sim->G, sim->particles[1], 0., a_dust, e_dust, inc_dust, Omega_dust, omega_dust, f_dust); 
    reb_add(sim, p2);

    /* REBOUNDx has a convenience function to calculate beta given the star's luminosity and the grain's physical radius, density and radiation pressure coefficient Q_pr (Burns et al. 1979). */
    
    double particle_radius = 1.e-5; /* in meters */
    double density = 1.e3;          /* kg/m3 = 1g/cc */
    double Q_pr = 1.;               /* Equals 1 in limit where particle radius >> wavelength of radiation */
    double L = 3.85e26;             /* Luminosity of the sun in Watts */

    beta = rebx_rad_calc_beta(rebx, params, particle_radius, density, Q_pr, L);
    rebx_set_param_double(&sim->particles[3], "beta", beta);    

    printf("Particle has beta = %f\n", rebx_get_param_double(&sim->particles[3], "beta"));

    reb_move_to_com(sim);

    reb_integrate(sim, tmax);
    rebx_free(rebx);                /* free memory allocated by REBOUNDx */
}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 1.e7)){
        const struct reb_particle* particles = sim->particles;
        const struct reb_particle saturn = particles[1];
        const double t = sim->t;
        const int N = sim->N;
        for (int i=2;i<N;i++){
            const struct reb_orbit orbit = reb_tools_particle_to_orbit(sim->G, sim->particles[i], saturn); /* calculate orbit of particles[i] around Saturn */
            printf("%e\t%e\t%e\n",t,orbit.a, orbit.e);
        }
    }
}
