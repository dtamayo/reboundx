/**
 * Self-Consistent Spin-Tidal and Dynamical equations of motion (Eggleton et. al 1998)
 *
 * In particular, this simulates the obliquity sculpting of the Kepler multis due to convergent migration
 * Result in based on Figure 3 in Millholland & Laughlin (2019), and reproduces Figure 3 in Lu et. al (2023)
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "tides_spin.c"

void heartbeat(struct reb_simulation* sim);

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    const double solar_mass = 1.;
    const double solar_rad = 0.00465;
    reb_add_fmt(sim, "m r", solar_mass, solar_rad);// Central object

    const double p1_mass = 5. * 3.0e-6; // in Earth masses * 1 Earth Mass / 1 Solar Mass
    const double p1_rad = 2.5 * 4.26e-5; // in Earth rad * 1 Earth rad / 1 AU
    reb_add_fmt(sim, "m a e r inc Omega pomega M", p1_mass, 0.17308688, 0.01, p1_rad, 0.5 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.)); // Planet 1

    const double p2_mass = 5. * 3.0e-6;
    const double p2_rad = 2.5 * 4.26e-5;
    reb_add_fmt(sim, "m a e r inc Omega pomega M", p2_mass, 0.23290608, 0.01, p2_rad, -0.431 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.)); // Planet 2
    reb_move_to_com(sim);
    sim->N_active = 3;
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->dt = 1e-3;

    // Add REBOUNDx Additional effects
    // First Spin

    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);
    // Sun
    const double solar_spin_period = 20 * 2 * M_PI / 365;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    const double solar_q = 100000000.;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", 0.07);
    //rebx_set_param_double(rebx, &sim->particles[0].ap, "sigma", 6303.);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "moi", 0.07 * solar_mass * solar_rad * solar_rad);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "spin_sx", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "spin_sy", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "spin_sz", solar_spin * 1.0);
    double solar_sigma = rebx_set_q(rebx, sim->G, &sim->particles[0], &sim->particles[1], solar_q);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sigma", solar_sigma);

    // P1
    const double spin_period_1 = 5. * 2. * M_PI / 365.; // 5 days in reb years
    const double spin_1 = (2. * M_PI) / spin_period_1;
    const double planet_q = 10000.;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.4);
    //rebx_set_param_double(rebx, &sim->particles[1].ap, "sigma", 1.75e15);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "moi", 0.25 * p1_mass * p1_rad * p1_rad);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "spin_sx", spin_1 * 0.0);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "spin_sy", spin_1 * -0.0261769);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "spin_sz", spin_1 * 0.99965732);
    double planet_sigma_1 = rebx_set_q(rebx, sim->G, &sim->particles[1], &sim->particles[0], planet_q);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sigma", planet_sigma_1);

    // P2
    double spin_period_2 = 3. * 2. * M_PI / 365.; // 3 days in reb years
    double spin_2 = (2. * M_PI) / spin_period_2;
    rebx_set_param_double(rebx, &sim->particles[2].ap, "k2", 0.4);
    // rebx_set_param_double(rebx, &sim->particles[2].ap, "sigma", 2.73e15);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "moi", 0.25 * p2_mass * p2_rad * p2_rad);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "spin_sx", spin_2 * 0.0);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "spin_sy", spin_2 * 0.0249736);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "spin_sz", spin_2 * 0.99968811);
    double planet_sigma_2 = rebx_set_q(rebx, sim->G, &sim->particles[2], &sim->particles[0], planet_q);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "sigma", planet_sigma_2);

    // And migration
    struct rebx_force* mo = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mo);

    // Set migration parameters
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -5e6 * 2 * M_PI);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", (-5e6 * 2 * M_PI) / 1.1);

    // Run simulation
    rebx_spin_initialize_ode(rebx, effect);

    FILE* f = fopen("test.txt","w");
    fprintf(f, "t,starx,stary,starz,starvx,starvy,starvz,star_sx, star_sy, star_sz, a1,i1,e1,s1x,s1y,s1z,mag1,pom1,Om1,f1,p1x,p1y,p1z,p1vx,p1vy,p1vz,a2,i2,e2,s2x,s2y,s2z,mag2,pom2,Om2,f2,p2x,p2y,p2z,p2vx,p2vy,p2vz\n");
    int cond = 0;
     for (int i=0; i<100000; i++){

         struct reb_particle* sun = &sim->particles[0];
         struct reb_particle* p1 = &sim->particles[1];
         struct reb_particle* p2 = &sim->particles[2];

  	 if (sim->t > 2e6 * 2 * M_PI && cond == 0){
  		printf("Migration Switching Off\n");
  		    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", INFINITY);
      		rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", INFINITY);
  	 	cond = 1;
  	 }

         double* star_sx = rebx_get_param(rebx, sun->ap, "spin_sx");
         double* star_sy = rebx_get_param(rebx, sun->ap, "spin_sy");
         double* star_sz = rebx_get_param(rebx, sun->ap, "spin_sz");

         double* sx1 = rebx_get_param(rebx, p1->ap, "spin_sx");
         double* sy1 = rebx_get_param(rebx, p1->ap, "spin_sy");
         double* sz1 = rebx_get_param(rebx, p1->ap, "spin_sz");

         double* sx2 = rebx_get_param(rebx, p2->ap, "spin_sx");
         double* sy2 = rebx_get_param(rebx, p2->ap, "spin_sy");
         double* sz2 = rebx_get_param(rebx, p2->ap, "spin_sz");

         struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, *p1, *sun);
         double a1 = o1.a;//vis_viva(r, &p1, &sun);
         double Om1 = o1.Omega;
         double i1 = o1.inc;
         double pom1 = o1.pomega;
         double f1 = o1.f;
         double e1 = o1.e;

         struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, *p2, *sun);
         double a2 = o2.a;//vis_viva(r, &p2, &sun);
         double Om2 = o2.Omega;
         double i2 = o2.inc;
         double pom2 = o2.pomega;
         double f2 = o2.f;
         double e2 = o2.e;

         struct reb_vec3d s1 = {*sx1, *sy1, *sz1};
         struct reb_vec3d s2 = {*sx2, *sy2, *sz2};

	 //printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", sim->t, a1, i1, e1, s1_inv.x, s1_inv.y, s1_inv.z, sqrt(s1_inv.x * s1_inv.x + s1_inv.y * s1_inv.y + s1_inv.z * s1_inv.z), pom1, Om1, p1->x, p1->y, p1->z, a2, i2, e2, s2_inv.x, s2_inv.y, s2_inv.z, sqrt(s2_inv.x * s2_inv.x + s2_inv.y * s2_inv.y + s2_inv.z * s2_inv.z), pom2, Om2, p2->x, p2->y, p2->z);


  //       struct reb_vec3d s1 = rebx_transform_inv_to_planet(i1, Om1, s1_inv);
  //       struct reb_vec3d s2 = rebx_transform_inv_to_planet(i2, Om2, s2_inv);

         // Interpret in the planet frame
         double mag1 = sqrt(s1.x * s1.x + s1.y * s1.y + s1.z * s1.z);
         double ob1 = acos(s1.z / mag1) * (180 / M_PI);
         double mag2 = sqrt(s2.x * s2.x + s2.y * s2.y + s2.z * s2.z);
         double ob2 = acos(s2.z / mag2) * (180 / M_PI);

         if (i % 100 == 0){
             printf("t=%f\t a1=%.6f\t a2=%.6f\t o1=%0.5f\t o2=%0.5f\n", sim->t / (2 * M_PI), a1, a2, ob1, ob2);
         }
         fprintf(f, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%0.10f,%0.10f,%0.10f,%0.10f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%0.10f,%0.10f,%0.10f,%0.10f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", sim->t / (2 * M_PI), sun->x, sun->y, sun->z, sun->vx, sun->vy, sun->vz, *star_sx, *star_sy, *star_sz, a1, i1, e1, s1.x, s1.y, s1.z, mag1, pom1, Om1, f1, p1->x,p1->y,p1->z, p1->vx, p1->vy, p1->vz, a2, i2, e2, s2.x, s2.y, s2.z, mag2, pom2, Om2, f2, p2->x,p2->y,p2->z, p2->vx, p2->vy, p2->vz);
         reb_integrate(sim,sim->t+(40 * 2 * M_PI));
     }
    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* r){
    if(reb_output_check(r, 25)){        // outputs to the screen
        //reb_output_timing(r, 1e4);
    }
}
