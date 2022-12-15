/**
 * Self-Consistent Spin-Tidal and Dynamical equations of motion (Eggleton et. al 1998)
 *
 * In particular, this simulates the pseudo-synchronization of a fiducial hot Jupiter
 * Definitely see the corresponding ipython example, as well as the documentation, for more explanations along the way of the various parameters and assumptions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "spin.c"

void heartbeat(struct reb_simulation* sim);
double tmax = 10000 * 2 * M_PI;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    const double solar_mass = 1.;
    const double solar_rad = 0.00465;
    reb_add_fmt(sim, "m r", solar_mass, solar_rad);// Central object

    const double p1_mass = 1. * 9.55e-4; // in Jupiter masses * 1 Jupiter Mass / 1 Solar Mass
    const double p1_rad = 1. * 4.676e-4; // in Jupiter rad * 1 jupiter rad / 1 AU
    const double p1_e = 0.01;

    reb_add_fmt(sim, "m a e r", p1_mass, 0.04072, p1_e, p1_rad); // Planet 1

    reb_move_to_com(sim);
    sim->N_active = 2;
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->dt = 1e-3;
    sim->heartbeat = heartbeat;

    // Add REBOUNDx Additional effects
    // First Spin
    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* effect = rebx_load_force(rebx, "spin");
    rebx_add_force(rebx, effect);
    // Star
    const double solar_spin_period = 27 * 2 * M_PI / 365;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    const double solar_q = 1000000.;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", 0.07);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "moi", 0.07 * solar_mass * solar_rad * solar_rad);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "spin_sx", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "spin_sy", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "spin_sz", solar_spin * 1.0);
    rebx_set_q(rebx, sim->G, &sim->particles[0], &sim->particles[1], solar_q);

    // Planet
    const double spin_period_1 = 0.5 * 2. * M_PI / 365.; // 0.5 days in reb years
    const double spin_1 = (2. * M_PI) / spin_period_1;
    const double planet_q = 10000.;
    const double theta_1 = 30. * (M_PI / 180.);
    const double phi_1 = 0 * (M_PI / 180);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.3);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "moi", 0.25 * p1_mass * p1_rad * p1_rad);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "spin_sx", spin_1 * sin(theta_1) * sin(phi_1));
    rebx_set_param_double(rebx, &sim->particles[1].ap, "spin_sy", spin_1 * sin(theta_1) * cos(phi_1));
    rebx_set_param_double(rebx, &sim->particles[1].ap, "spin_sz", spin_1 * cos(theta_1));
    rebx_set_q(rebx, sim->G, &sim->particles[1], &sim->particles[0], planet_q);

    rebx_align_simulation(rebx);
    rebx_spin_initialize_ode(rebx, effect);
/*
    FILE* f = fopen("12_12_hj_spindown_reb_updated.txt","w");
    fprintf(f,"t,a1,i1,e1,s1x,s1y,s1z,mag1,pom1,Om1\n");
    //printf("t,starx,stary,starz,starvx,starvy,starvz,star_sx,star_sy,star_sz,a1,i1,e1,s1x,s1y,s1z,mag1,pom1,Om1,f1,p1x,p1y,p1z,p1vx,p1vy,p1vz\n");
    /*
    struct reb_particle* star = &sim->particles[0];
    struct reb_particle* p = &sim->particles[1];

    struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, *p, *star);
    double a = o.a;//vis_viva(r, &p1, &sun);
    double Om = o.Omega;
    double i = o.inc;
    double pom = o.pomega;
    double e = o.e;

     for (int i=0; i<100000; i++){

         struct reb_particle* sun = &sim->particles[0];
         struct reb_particle* p1 = &sim->particles[1];

         double* sx1 = rebx_get_param(rebx, p1->ap, "spin_sx");
         double* sy1 = rebx_get_param(rebx, p1->ap, "spin_sy");
         double* sz1 = rebx_get_param(rebx, p1->ap, "spin_sz");

         struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, *p1, *sun);
         double a1 = o1.a;//vis_viva(r, &p1, &sun);
         double Om1 = o1.Omega;
         double i1 = o1.inc;
         double pom1 = o1.pomega;
         double e1 = o1.e;

         struct reb_vec3d s1 = {*sx1, *sy1, *sz1};

         // Interpret in the planet frame
         double mag1 = sqrt(s1.x * s1.x + s1.y * s1.y + s1.z * s1.z);
         double ob1 = acos(s1.z / mag1) * (180 / M_PI);

	       if (i % 1000 == 0){
             printf("t=%f\t a1=%.6f\t o1=%0.5f, mag1=%0.5f\t", sim->t / (2 * M_PI), a1, ob1, mag1);
         }
        fprintf(f, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", sim->t / (2 * M_PI), a1, i1, e1, s1.x, s1.y, s1.z, mag1, pom1, Om1);
        reb_integrate(sim, sim->t+(1 * 2 * M_PI));
     }
     */
    system("rm -v output.txt");
    reb_integrate(sim, tmax);
    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Output spin and orbital information to file
    if(reb_output_check(sim, 10)){        // outputs every 10 REBOUND years
      struct rebx_extras* const rebx = sim->extras;
      FILE* of = fopen("output.txt", "a");
      if (of==NULL){
          reb_error(sim, "Can not open file.");
          return;
      }

      struct reb_particle* star = &sim->particles[0];
      struct reb_particle* p = &sim->particles[1];

      double* sx = rebx_get_param(rebx, p->ap, "spin_sx");
      double* sy = rebx_get_param(rebx, p->ap, "spin_sy");
      double* sz = rebx_get_param(rebx, p->ap, "spin_sz");
      double mag = sqrt(*sx * *sx + *sy * *sy + *sz * *sz);
      double ob = acos(*sz / mag) * (180 / M_PI);

      struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, *p, *star);
      double a = orb.a;
      double Om = orb.Omega;
      double inc = orb.inc;
      double pom = orb.pomega;
      double e = orb.e;

      fprintf(of, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", sim->t, a, inc, e, mag, pom, Om, *sx, *sy, *sz);

      fclose(of);
    }

    if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
        reb_output_timing(sim, tmax);
    }
}
