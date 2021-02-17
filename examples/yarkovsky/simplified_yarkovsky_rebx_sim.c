//
//  simplified_yarkovsky_rebx_sim.c
//  
//
//  Created by Noah Ferich on 2/13/21.
//

#include "simplified_yarkovsky_rebx_sim.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]) {

struct reb_simulation* sim = reb_create_simulation(); //creates simulation

sim->G = 4*M_PI*M_PI;  // use units of AU, yr and solar masses
sim->dt = .05;         //timestep for simulation in yrs
sim->integrator = REB_INTEGRATOR_WHFAST; //integrator for sim

//following adds star with mass of Sun to sim
struct reb_particle star = {0};
star.m = 1.;
reb_add(sim, star);

//following variables are the orbital elements of only asteroid in sim
double m = 0;
double a = .5;
double e = 0;
double inc = 0;
double Omega = 0;
double omega = 0;
double f = 0;

//adds asteroid to the sim
struct reb_particle asteroid_1 = reb_tools_orbit_to_particle(sim->G, star, m, a, e, inc, Omega, omega, f);
reb_add(sim,asteroid_1);

double a_2 = .75;
struct reb_particle asteroid_2 = reb_tools_orbit_to_particle(sim->G, star, m, a_2, e, inc, Omega, omega, f);
reb_add(sim,asteroid_2);
        
double a_3 = 1.0;
    
struct reb_particle asteroid_3 = reb_tools_orbit_to_particle(sim->G, star, m, a_3, e, inc, Omega, omega, f);
reb_add(sim,asteroid_3);

    
struct reb_particle* const particles = sim->particles; //pointer for the particles in the sim

struct rebx_extras* rebx = rebx_attach(sim);
struct rebx_force* yark = rebx_load_force(rebx, "yark_max_out");
rebx_set_param_double(rebx, &sim->particles[1].ap, "body_density", 3000);
rebx_set_param_double(rebx, &yark->ap, "lstar", 3.828e26);
    particles[1].r = 1000;
    
rebx_set_param_double(rebx, &sim->particles[3].ap, "body_density", 3000);
//rebx_set_param_double(rebx, &sim->particles[3].ap, "lstar", 3.828e26);
   particles[3].r = 1000;

rebx_add_force(rebx, yark);

    double tmax = 50000;
    
reb_integrate(sim, tmax); //integrates system for tmax years
    
struct reb_orbit o= reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]); //o gives orbital parameters for asteroid after sim
struct reb_orbit p= reb_tools_particle_to_orbit(sim->G, sim->particles[2], sim->particles[0]); //o gives orbital parameters for asteroid after sim
struct reb_orbit q= reb_tools_particle_to_orbit(sim->G, sim->particles[3], sim->particles[0]); //o gives orbital parameters for asteroid after sim
    
double final_a = o.a; //final semi-major axis of asteroid after sim
    
    double final_a_2 = p.a;
    
    double final_a_3 = q.a;
    
printf("CHANGE IN SEMI-MAJOR AXIS: %1.30f\n", (final_a-a)); //prints difference between the intitial and final semi-major axes of asteroid

printf("CHANGE IN SEMI-MAJOR AXIS: %1.30f\n", (final_a_2-a_2));
    
printf("CHANGE IN SEMI-MAJOR AXIS: %1.30f\n", (final_a_3-a_3));
    
}
