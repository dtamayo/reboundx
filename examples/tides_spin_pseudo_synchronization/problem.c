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
#include "tides_spin.c"

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

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
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
    double solar_sigma = rebx_set_q(rebx, sim->G, &sim->particles[0], &sim->particles[1], solar_q);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sigma", solar_sigma);
    //rebx_set_q(rebx, sim->G, &sim->particles[0], &sim->particles[1], solar_q);

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
    double planet_sigma = rebx_set_q(rebx, sim->G, &sim->particles[1], &sim->particles[0], planet_q);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sigma", planet_sigma);
    // rebx_set_q(rebx, sim->G, &sim->particles[1], &sim->particles[0], planet_q);

    rebx_align_simulation(rebx);
    rebx_spin_initialize_ode(rebx, effect);

    system("rm -v output.txt"); // remove previous output file
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
      // double ob = acos(*sz / mag) * (180 / M_PI);

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
