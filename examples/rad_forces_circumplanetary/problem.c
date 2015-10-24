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

double tmax = 1e11;

int main(int argc, char* argv[]){
	struct reb_simulation* sim = reb_create_simulation();
	// Setup constants
	sim->integrator		= REB_INTEGRATOR_IAS15;
	sim->G				= 6.674e-11;	// Use SI units
	sim->dt 			= 1e4;			// Initial timestep in sec (~1 hr.  IAS15 will adjust adaptively)
	sim->N_active		= 2;			// Only the sun and the planet affect other particles gravitationally (simulation will be very slow without this!)
	//sim->heartbeat	 	= heartbeat;
	sim->usleep		= 5000;				// Slow down integration (for visualization only)

	struct rebx_extras* rebx = rebx_init(sim); // Always add rebx to simulation first
	double c = 3.e8;					// speed of light in SI units
	double L = 3.85e26;					// Luminosity of the sun in SI units (W)
	rebx_add_radiation_forces(rebx, c, L); // add radiation forces
	
	// sun
	struct reb_particle sun = {0};
	sun.m  = 1.99e30;					// mass of Sun in kg
	reb_add(sim, sun);

	// Saturn (simulation set up in Saturn's orbital plane, i.e., inc=0, so only need 4 orbital elements
	double mass_sat = 5.68e26;			// Mass of Saturn
	double a_sat = 1.43e12;				// Semimajor axis of Saturn in m
	double e_sat = 0.056;				// Eccentricity of Saturn
	double pomega_sat = 0.;				// Angle from x axis to pericenter
	double f_sat = 0.;					// True anomaly of Saturn
	struct reb_particle saturn = reb_tools_orbit2d_to_particle(sim->G, sun, mass_sat, a_sat, e_sat, pomega_sat, f_sat);

	reb_add(sim, saturn);

	/* Dust particles
	 Here we imagine particles launched from Saturn's irregular Satellite Phoebe
	 Since the escape velocity from the moon << Phoebe's orbital velocity, we
	 can assume particles initially inherit the moon's position and velocity (e.g. Tamayo et al. 2011)
	 We therefore initialize the dust particle orbits here with Phoebe's current osculating 
	 orbital elements (though we set most of the angles to 0).  The dynamics will be different 
	 depending on these angles as well as the Sun's initial position.*/
		
	/* In order for a particle to feel radiation forces, we have to set the radiation pressure 
	 * parameter Q_pr, as well as the particle's mass and radius.  We do this in two ways below:
	 * in both cases we use the particle's bulk density, and then we use either the particle's 
	 * radius directly, or calculate it given a value of $\beta, the ratio of the radiation
	 * pressure force to the gravitional force from the star (Burns et al. 1979). */
		
	double a_dust = 1.30e10;// semimajor axis of satellite Phoebe, in m
	double e_dust = 0.16;	// eccentricity of Phoebe
	double inc_dust = 175.*M_PI/180.;	// inclination of Phoebe to Saturn's orbital plane
	double Omega_dust = 0.;	// longitude of ascending node
	double omega_dust = 0.; // argument of pericenter
	double f_dust = 0.;		// true anomaly

	// we first set up the orbit
	struct reb_particle p = reb_tools_orbit_to_particle(sim->G, sim->particles[1], 0., a_dust, e_dust, inc_dust, Omega_dust, omega_dust, f_dust); 
	// we pass 0 mass for the particle in the 3rd parameter here for simplicity since it won't affect the orbit.

	// now we set the physical parameters
	
	double Q_pr = 1.;				// Radiation pressure coefficient. Equals 1 in limit where particle radius >> wavelength of radiation
	rebx_set_Q_pr(rebx, &p, Q_pr); 	// Only particles with Q_pr set will feel radiation forces.
	printf("%p\n", &p.ap);
	double density_dust = 1.e3;		// kg/m^3 = 1g/cc
	p.r = 1.e-5;					// dust grain radius in m 
	p.m = rebx_rad_calc_mass(density_dust, p.r);	// assumes spherical grains
	reb_add(sim, p); 
	printf("%p\t%p\n", &p.ap, &sim->particles[2].ap);
	rebx_set_Q_pr(rebx, &sim->particles[2], Q_pr);	
	
	// Now add a 2nd particle of different size and density on same orbit, using a beta value
	struct reb_particle p2 = reb_tools_orbit_to_particle(sim->G, sim->particles[1], 0., a_dust, e_dust, inc_dust, Omega_dust, omega_dust, f_dust); 

	//rebx_set_Q_pr(rebx, &p2, Q_pr);
	density_dust = 3.e3;							
	double beta = 1.e-3;
	p2.r = rebx_rad_calc_particle_radius(rebx, beta, density_dust, Q_pr);
	p2.m = rebx_rad_calc_mass(density_dust, p.r);	// assumes spherical grains
	reb_add(sim, p2); 

	reb_move_to_com(sim);

	system("rm -v a.txt");	

	printf("m\tr\tbeta\n");
	printf("%e\t%e\t%e\n", sim->particles[2].m, sim->particles[2].r, rebx_rad_calc_beta(rebx, &sim->particles[2]));
	printf("%e\t%e\t%e\n", rebx_rad_calc_mass(1.e3, sim->particles[2].r), sim->particles[2].r, rebx_rad_calc_beta(rebx, &sim->particles[2]));

	reb_integrate(sim, tmax);
}

void heartbeat(struct reb_simulation* sim){
	if(reb_output_check(sim, 1.e8)){
		reb_output_timing(sim, tmax);
	}
	if(reb_output_check(sim, 1.e8)){ // output every year
		FILE* f = fopen("a.txt","a");
		const struct reb_particle* particles = sim->particles;
		const struct reb_particle saturn = particles[1];
		const double t = sim->t;
		const int N = sim->N;
		for (int i=2;i<N;i++){
			const struct reb_orbit orbit = reb_tools_particle_to_orbit(sim->G, sim->particles[i], saturn); // calculate orbit of particles[i] around Saturn
			fprintf(f,"%e\t%e\t%e\n",t,orbit.a, orbit.e);
		}
		fclose(f);
	}
}
