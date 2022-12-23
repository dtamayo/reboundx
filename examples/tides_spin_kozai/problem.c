/**
 * Kozai cycles
 *
 * This example uses the IAS15 integrator to simulate
 * a Lidov Kozai cycle of a planet perturbed by a distant star.
 * The integrator automatically adjusts the timestep so that
 * even very high eccentricity encounters are resolved with high
 * accuracy.
 *
 * This is the same Kozai example implemented in base REBOUND, modified to include the tidal and spin effects from bodies with structure.
 * Please refer to that example for the system details
 * Also, see the ipython examples prefixed TidesSpin for in-depth exploration of the parameters that can be set in this simulation.
 */
 #include <stdio.h>
 #include <stdlib.h>
 #include <unistd.h>
 #include <math.h>
 #include "rebound.h"
 #include "reboundx.h"
 #include "tides_spin.c"

void heartbeat(struct reb_simulation* r);
double tmax = 4e3; // kept short to run quickly. set to at least 4x longer to see full evolution

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Initial conditions
    // Setup constants
    sim->dt             = M_PI*1e-2;     // initial timestep
    sim->integrator        = REB_INTEGRATOR_IAS15; // IAS15 is used for its adaptive timestep:
                                                   // in a Kozai cycle the planet experiences close encounters during the high-eccentricity epochs.
                                                   // A fixed-time integrator (for example, WHFast) would need to apply the worst-case timestep to the whole simulation
    sim->heartbeat        = heartbeat;

    // Initial conditions
    struct reb_particle star = {0};
    star.m  = 1;
    star.r = 0.00465;
    reb_add(sim, star);

    struct reb_particle planet = {0};
    planet.m  = 0.05 * 9.55e-4; // A Neptune-like planet
    planet.r = 0.3 * 4.676e-4;
    double e_planet = 0;
    planet.x  = 1. - e_planet;
    planet.vy = sqrt((1. + e_planet) / (1. - e_planet));
    reb_add(sim, planet);

    // The perturber - treated as a point particle
    struct reb_particle perturber = {0};
    perturber.x  = 10;
    double inc_perturber = 89.9;
    perturber.m  = 1;
    perturber.vy = cos(inc_perturber/180.*M_PI)*sqrt((star.m+perturber.m)/perturber.x);
    perturber.vz = sin(inc_perturber/180.*M_PI)*sqrt((star.m+perturber.m)/perturber.x);
    reb_add(sim, perturber);

    // Add REBOUNDx effects
    // First tides_spin
    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    // Sun
    const double solar_spin_period = 27 * 2. * M_PI / 365.;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    const double solar_k2 = 0.1;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "moi", 0.07 * star.m * star.r * star.r);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sx", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sy", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sz", solar_spin * 1.0);
    double solar_sigma = rebx_tides_calc_sigma_from_Q(rebx, &sim->particles[0], &sim->particles[1], 3e6);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sigma", solar_sigma);

    // P1
    const double spin_period_p = 1. * 2. * M_PI / 365.; // days to reb years
    const double spin_p = (2. * M_PI) / spin_period_p;
    const double planet_k2 = 0.4;
    const double planet_q = 1e4;
    const double theta_1 = 0. * M_PI / 180.;
    const double phi_1 = 0. * M_PI / 180;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", planet_k2);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "moi", 0.25 * planet.m * planet.r * planet.r);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sx", spin_p * sin(theta_1) * sin(phi_1));
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sy", spin_p * sin(theta_1) * cos(phi_1));
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sz", spin_p * cos(theta_1));
    double planet_sigma = rebx_tides_calc_sigma_from_Q(rebx, &sim->particles[1], &sim->particles[0], planet_q);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sigma", planet_sigma);


    // add GR precession:
    struct rebx_force* gr = rebx_load_force(rebx, "gr_full");
    rebx_add_force(rebx, gr);
    rebx_set_param_double(rebx, &gr->ap, "c", 10065.32); // in default units

    reb_move_to_com(sim);
    rebx_align_simulation(rebx);
    rebx_spin_initialize_ode(rebx, effect);

    system("rm -v orbits.txt");        // delete previous output file
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

      struct reb_particle* sun = &sim->particles[0];
      struct reb_particle* p1 = &sim->particles[1];
      struct reb_particle* pert = &sim->particles[2];

      double* sx_sun = rebx_get_param(rebx, sun->ap, "sx");
      double* sy_sun = rebx_get_param(rebx, sun->ap, "sy");
      double* sz_sun = rebx_get_param(rebx, sun->ap, "sz");
      double mag_sun = sqrt((*sx_sun) * (*sx_sun) + (*sy_sun) * (*sy_sun) + (*sz_sun) * (*sz_sun));

      double* sx_p = rebx_get_param(rebx, p1->ap, "sx");
      double* sy_p = rebx_get_param(rebx, p1->ap, "sy");
      double* sz_p = rebx_get_param(rebx, p1->ap, "sz");
      double mag_p = sqrt((*sx_p) * (*sx_p) + (*sy_p) * (*sy_p) + (*sz_p) * (*sz_p));

      struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, *p1, *sun);
      double a1 = o1.a;
      double Om1 = o1.Omega;
      double i1 = o1.inc;
      double pom1 = o1.pomega;
      double e1 = o1.e;

      struct reb_particle com = reb_get_com_of_pair(sim->particles[0],sim->particles[1]);
      struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, *pert, com);
      double a2 = o2.a;
      double Om2 = o2.Omega;
      double i2 = o2.inc;
      double pom2 = o2.pomega;
      double e2 = o2.e;

      fprintf(of, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", sim->t, sx_sun, sy_sun, sz_sun, mag_sun, sx_p, sy_p, sz_p, mag_p, a1, Om1, i1, pom1, e1, a2, Om2, i2, pom2, e2); // print spins and orbits

      fclose(of);
    }

    if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
        reb_output_timing(sim, tmax);
    }
}
