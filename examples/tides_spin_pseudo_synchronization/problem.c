/**
 * Self-Consistent Spin-Tidal and Dynamical equations of motion (Eggleton et. al 1998)
 *
 * This example shows how to add quadrupole and tidal distortions to bodies with structure, letting us consistently track the spin-axis and dynamical evolution of the system.
 * In particular, this simulates the pseudo-synchronization of a fiducial hot Jupiter.
 * The hot Jupiter is initialized with a slight eccentricity, nontrivial obliquity and fast rotation.
 * Under the influence of tidal dissipation, we see the following rapidly occur: circularization of the orbit, obliquity damping down to 0, and rotation settling to the pseudo-synchronous value predicted by Hut (1981)
 * Definitely see the corresponding ipython example, as well as the documentation, for more in-depth explanations regarding the various parameters that may be set in this simulation.
 * Also see Lu et. al (in review), Eggleton et. al (1998), Hut (1981).
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

    // Star
    const double solar_mass = 1.;
    const double solar_rad = 0.00465; // Bodies with structure require radius! This is one solar radius
    reb_add_fmt(sim, "m r", solar_mass, solar_rad);// Central object

    // Fiducual hot Jupiter
    const double p1_mass = 1. * 9.55e-4; // in Jupiter masses * 1 Jupiter Mass / 1 Solar Mass
    const double p1_rad = 1. * 4.676e-4; // in Jupiter rad * 1 jupiter rad / 1 AU
    const double p1_e = 0.01;
    reb_add_fmt(sim, "m a e r", p1_mass, 0.04072, p1_e, p1_rad); // Planet 1

    sim->N_active = 2;
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->dt = 1e-3;
    sim->heartbeat = heartbeat;

    // Add tides_spin as a REBOUNDx additional effect
    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);
    // Star
    // The following parameters are mandatory to set for the body to have structure:
    // The Love number (k2)
    // The three components of spin frequency (sx, sy, sz)
    const double solar_k2 = 0.07;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    const double solar_spin_period = 27 * 2 * M_PI / 365; // 27 Days in REBOUND time units
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sx", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sy", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sz", solar_spin * 1.0);

    // While not technically necessary, the fully dimensional moment of inertia (moi) should also be set
    // This is required to evolve the spin axis! Without setting this value, the spin axis will remain stationary.
    rebx_set_param_double(rebx, &sim->particles[0].ap, "moi", 0.07 * solar_mass * solar_rad * solar_rad);

    // Finally, the last parameter is the tidal dissipation parameter (sigma)
    // You can either set sigma directly, or use two helper functions that relate sigma to more intuitive quantities:
    // 'rebx_tides_calc_sigma_from_Q' uses the Tidal Quality Factor (Q)
    // 'rebx_tides_calc_sigma_from_tau' uses the constant time lag (tau)
    // See Lu et. al (in review) for the specific relations used to relate these quantities.
    const double solar_Q = 1000000.;

    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    const double solar_tau = 1 / (2 * solar_Q * orb.n);

    // You can set sigma directly whenever you want, but the two helper functions depend on the previously set values above
    // For this reason, if you choose to use the helper functions, USE THEM LAST, after the other tides_spin parameters have been set!
    // Here are examples of all three ways to set sigma. All three are equivalent!

    // rebx_tides_calc_sigma_from_Q takes, as input:
    // 1) The rebx extras object
    // 2) The simulation's value of G
    // 3) The particle you wish to set sigma for
    // 4) A tidal perturber
    // 5) The Q value you will use
    double solar_sigma = rebx_tides_calc_sigma_from_Q(rebx, sim->G, &sim->particles[0], &sim->particles[1], solar_Q);

    // rebx_tides_calc_sigma_from_tau takes, as input:
    // 1) The rebx extras object
    // 2) The simulation's value of G
    // 3) The particle you wish to set sigma for
    // 4) The Q value you will use
    // double solar_sigma = rebx_tides_calc_sigma_from_tau(rebx, sim->G, &sim->particles[0], solar_tau);

    //double solar_sigma = 4 * sim->G * solar_tau / (3 * solar_rad * solar_rad * solar_rad * solar_rad * solar_rad * solar_k2);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sigma", solar_sigma);

    // Planet - all of the above applies here too
    const double spin_period_1 = 0.5 * 2. * M_PI / 365.; // 0.5 days in reb years
    const double spin_1 = (2. * M_PI) / spin_period_1;
    const double planet_Q = 10000.;
    const double theta_1 = 30. * (M_PI / 180.);
    const double phi_1 = 0 * (M_PI / 180);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.3);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "moi", 0.25 * p1_mass * p1_rad * p1_rad);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sx", spin_1 * sin(theta_1) * sin(phi_1));
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sy", spin_1 * sin(theta_1) * cos(phi_1));
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sz", spin_1 * cos(theta_1));
    double planet_sigma = rebx_tides_calc_sigma_from_Q(rebx, sim->G, &sim->particles[1], &sim->particles[0], planet_Q);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sigma", planet_sigma);

    reb_move_to_com(sim);
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

      double* sx = rebx_get_param(rebx, p->ap, "sx");
      double* sy = rebx_get_param(rebx, p->ap, "sy");
      double* sz = rebx_get_param(rebx, p->ap, "sz");
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
