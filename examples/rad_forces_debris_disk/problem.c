/**
 * Radiation forces on a debris disk
 * 
 * This example shows how to integrate a debris disk around a Sun-like star, with
 * dust particles under the action of radiation forces using WHFast. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* r);

double tmax = 3e12;	// in sec, ~ 10^5 yrs

int main(int argc, char* argv[]){
	struct reb_simulation* sim = reb_create_simulation();
	// Setup constants
	double AU = 1.5e11;					// in m
	sim->integrator		= REB_INTEGRATOR_WHFAST;
	sim->G				= 6.674e-11;	// Use SI units
	sim->dt 			= 1e8;			// At ~100 AU, orbital periods are ~1000 yrs, so here we use ~1% of that, in sec
	sim->N_active		= 1;			// The dust particles do not interact with one another gravitationally
	sim->heartbeat	 	= heartbeat;

	struct rebx_extras* rebx = rebx_init(sim); 
	double c = 3.e8;					// speed of light in SI units
	double L = 3.85e26;					// Luminosity of the sun in SI units (W)
	rebx_add_radiation_forces(rebx, c, L); // add radiation forces
	
	// sun
	struct reb_particle sun = {0};
	sun.m  = 1.99e30;					// mass of Sun in kg
	reb_add(sim, sun);

	/* Dust particles
	 * We idealize a perfectly coplanar debris disk with particles that have semimajor axes between 100 and 120 AU.
	 * We initialize particles with 0.01 eccentricity, random pericenters and azimuths, and uniformly distributed
	 * semimajor axes in the above range.  We take all particles to have a beta parameter (ratio of radiation 
	 * force to gravitational force from the star) of 0.1, and a density of 1g/cc.
	 *
	 * We also have to set each particle's radidation pressure coefficient Q_pr (Burns et al. 1979).  Only particles
	 * with Q_pr set will feel radiation forces.  For particles with radii >> the stellar radiation wavelength, Q_pr = 1.*/

	double amin = 100.*AU;
	double awidth = 20.*AU;
	double e = 0.01;
	double Ndust = 1000;				// Number of dust particles

	int seed = 3;					// random number generator seed
	srand(seed);
	double a, pomega, f;
	struct reb_particle p;
	double beta = 0.1;
	double dust_density = 1000.;	// kg/m^3 = 1g/cc
	double Q_pr = 1.;				
	for(int i=1; i<=Ndust; i++){
		/* first we set up the orbit.  For coplanar orbits, we can use reb_tools_orbit2d_to_particle to initialize
		 a reb_particle from a set of orbital elements.  For simplicity, we pass it a dust particle mass of 0
		 since it won't matter for the orbit.*/
		a = amin + awidth*(double)rand() / (double)RAND_MAX;
		pomega = 2*M_PI*(double)rand() / (double)RAND_MAX;
		f = 2*M_PI*(double)rand() / (double)RAND_MAX;
			
		p = reb_tools_orbit2d_to_particle(sim->G, sim->particles[0], 0., a, e, pomega, f);
		
		/* Now we set the physical parameters.  For radiation forces, we have to set the particle's mass and radius.
		 Here we choose to parametrize things in terms of a bulk density and beta parameter, and use REBOUNDx 
		 convenience functions to calculate the mass and radius.*/
		
		p.r = rebx_rad_calc_particle_radius(rebx, beta, dust_density, Q_pr); 	// convenience function for getting size from beta
		p.m = rebx_rad_calc_mass(dust_density, p.r);							// assumes spherical grains
		reb_add(sim, p); 
		
		/*We can only call rebx parameter setters on particles in the sim->particles array, so we do this after adding the particle to sim.*/
		rebx_set_Q_pr(&sim->particles[i], Q_pr); 	// Only particles with Q_pr set will feel radiation forces.
	}

	reb_move_to_com(sim);

	reb_integrate(sim, tmax);
	rebx_free(rebx);								// free memory allocated by REBOUNDx
}

void heartbeat(struct reb_simulation* sim){
	if(reb_output_check(sim, 1.e8)){
		reb_output_timing(sim, tmax);
	}
	/* You could also write output to a file here.*/
}
