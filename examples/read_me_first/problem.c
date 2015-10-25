/**
 * General considerations for simulations using REBOUNDx
 *
 * This example goes through some of the general requirements
 * for all REBOUNDx simulations, some general recipes that
 * may be useful, and possible surprises to watch out for.
 *
 * For techniques general to REBOUND, see the rebound/examples
 * directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]){
	struct reb_simulation* sim = reb_create_simulation();

	/**************************
	 * Setting up the simulation
	 ***************************/

	/*Everything in the simulation has to be entered in the same set of units.  This can be tricky
	 * when adding additional effects, since the scales can be so different.  For example if the orbits
	 * are in AU, then dust particles need their radii also specified in AU!*/
	sim->G				= 6.674e-11;	// Here we choose to use SI units

	/*If you have a lot of low mass particles (planetesimals, dust grains etc.), you might want to have
	 * these interact with the higher-mass planets gravitationally, but not with one another, to speed
	 * up the code.  To do this you have to add the high-mass particles to the simulation FIRST, and
	 * you then set the N_active variable to the number of high-mass particles.*/
	sim->N_active		= 1;			// Set so none of the circumstellar bodies interact with each other gravitationally

	/*Before we call any REBOUNDx functions, we always have to add REBOUNDx to our simulation like this:*/
	struct rebx_extras* rebx = rebx_init(sim);
	
	/*Now we can add different effects to our simulation.  Some take different constants, which we have to
	 * pass since REBOUNDx doesn't know what units we want to use in our simulation. All functions to add
	 * effects are named in the format rebx_add_effect(rebx, constants)*/
	double c = 3.e8;						// speed of light in SI units
	double L = 3.85e26;						// Luminosity of the sun in SI units (W)
	rebx_add_radiation_forces(rebx, c, L); 	// add radiation forces
	rebx_add_gr(rebx, c);					// add post-newtonian corrections from general relativity
	rebx_add_modify_orbits_direct(rebx);	// add orbit modifications (doesn't take any parameters)

	/*At any time, we can change any of the parameters specific to a REBOUNDx modification by accessing
	 * the modification structure directly.  For a given effect and param, one accesses it through
	 * rebx->effect.param:*/

	rebx->radiation_forces.L = 8.e26;				// Change the luminosity of the central body
	printf("L = %e\n", rebx->radiation_forces.L);	// prints 8.e26

	/******************
	 * Adding particles
	 ******************/
	
	/* When you declare a reb_particle structure, all the contents of its variables are undefined.  To avoid
	 * remembering all the fields and manually setting them, the easiest way is to always initialize at
	 * declaration with {0}.  This sets all fields to 0, including the pointers (i.e. NULL).  If you do not
	 * set the pointers to NULL like this (or manually) REBOUNDx will likely crash.*/
	struct reb_particle sun = {0}; 		// Since all fields are 0, this puts Sun at rest at the origin
	sun.m  = 1.99e30;					// mass of Sun in kg
	reb_add(sim, sun);

	/*To keep the focus on REBOUNDx, we don't initialize the orbits of any of the particles below*/

	struct reb_particle p1 = {0};
	/*************************************
	 * Adding particle-specific parameters
	 *************************************/

	/* Different REBOUNDx modifications have particle-specific parameters.  These can be set by functions
	 * always of the form rebx_set_param(&particle, value) (note one must pass a pointer to the particle).
	 * There is also a corresponding getter rebx_get_param(&particle).
	 * IMPORTANT:  You can only call these functions with particles in the sim->particles array, i.e.
	 * you have to set up a reb_particle, add it to the simulation, and THEN add REBOUNDx parameters
	 * using a pointer to sim->particles[i].  rebx_set_tau_a(p1, -1.e7); would give an error.*/

	reb_add(sim, p1);
	rebx_set_tau_a(&sim->particles[1], -1.e7);	// add a semimajor axis damping timescale of 1.e7 sec.
	printf("p1's tau_a = %e\n", rebx_get_tau_a(&sim->particles[1]));	// prints -1.e7

	/* If you are adding lots of particles with lots of REBOUNDx parameters that all need to be the same, 
	 * you might want to clone lots of particles.  So what happens if we copy a particle?*/

	struct reb_particle p2 = sim->particles[1];
	reb_add(sim, p2);

	printf("p2's tau_a = %e\n", rebx_get_tau_a(&sim->particles[2]));	// This also prints -1.e7, so it seems to have worked, but

	rebx_set_tau_a(&sim->particles[2], -3.e7);	
	printf("p2's tau_a = %e\n", rebx_get_tau_a(&sim->particles[2]));	// Prints -3.e7 as expected
	printf("p1's tau_a = %e\n", rebx_get_tau_a(&sim->particles[1])); // Also prints -3.e7!

	/* This surprise happens because in C, assigning structures to one another as we did above makes a shallow copy. This
	 * means that all the fields in the structure get their own independent copy.  The problem is that if there are any 
	 * fields that are pointers, just the pointer gets copied, so now the two copies point at a single spot in memory.
	 * The particle REBOUNDx parameters are stored in the particle.ap pointer so that's why any changes in the REBOUNDx 
	 * parameters for particles[2] get reflected in particles[1] and vice-versa.  Whether that's a problem depends on your
	 * application, but it's important to be aware of this potential pitfall.*/

	/* Now you would integrate, etc. as normal */
}
