/**
 * Velocity dependent drag force
 *
 * This is a very simple example on how to implement a velocity 
 * dependent drag force. The example uses the IAS15 integrator, which 
 * is ideally suited to handle non-conservative forces.
 * No gravitational forces or collisions are present.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* const r);

double tmax = 1.e6;

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->dt 		= 0.012;		// initial timestep.
	r->integrator	= REB_INTEGRATOR_WHFAST;
	//r->integrator	= REB_INTEGRATOR_IAS15;

	struct reb_particle p = {0}; 
	p.m  	= 1.;	
	reb_add(r, p); 

	struct reb_particle p1 = reb_tools_orbit2d_to_particle(r->G, p,  1.e-8, 1.0, 0.4, 0., 0.);	
	reb_add(r,p1);

	struct rebx_extras* rebx = rebx_init(r);

	rebx_add_gr(r,10000.);
	rebx->gr.c /=100.; // enhance precession

	reb_move_to_com(r);

	reb_integrate(r, tmax);
}

void heartbeat(struct reb_simulation* const r){
	// Output some information to the screen every 100th timestep
	if(reb_output_check(r, 100.*r->dt)){
		//struct reb_orbit o1 = reb_tools_p2orbit(r->G, r->particles[1], r->particles[0]);
		//struct reb_orbit o2 = reb_tools_p2orbit(r->G, r->particles[2], r->particles[0]);
		//printf("%f\t%f\t%f\n", r->t, o1.a, o2.a);
		//reb_output_timing(r, tmax);
	}
}
