/**
 * Kozai cycles
 *
 * This example uses the IAS15 integrator to simulate
 * a Lidov Kozai cycle of a planet perturbed by a distant star.
 * The integrator automatically adjusts the timestep so that
 * even very high eccentricity encounters are resolved with high
 * accuracy.
 *
 * This is the same Kozai example implemented in Lu et. al (2023)
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
double tmax = 1e5; // kept short to run quickly.
                   // set to 3e5 for a full cycle
                   // or 7e6 * 2 * M_PI to reproduce the paper plot

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Initial conditions
    // Setup constants
    sim->dt                 = M_PI*1e-1;     // initial timestep
    sim->integrator         = REB_INTEGRATOR_IAS15; // IAS15 is used for its adaptive timestep:
                                                   // in a Kozai cycle the planet experiences close encounters during the high-eccentricity epochs.
                                                   // A fixed-time integrator (for example, WHFast) would need to apply the worst-case timestep to the whole simulation
    sim->heartbeat          = heartbeat;

    // Initial conditions
    struct reb_particle star = {0};
    star.m  = 1;
    star.r = 0.00465;
    reb_add(sim, star);

    // struct reb_particle planet = {0};
    double planet_m  = 0.054 * 9.55e-4; // A Jupiter-like planet
    double planet_r = 0.3 * 4.676e-4;
    double planet_a = 2.;
    double planet_e = 0.001;
    reb_add_fmt(sim, "m r a e", planet_m, planet_r, planet_a, planet_e);

    // The perturber - treated as a point particle
    double perturber_m  = 1;
    double perturber_a = 50.;
    double perturber_e = 0.7 * M_PI / 180.;
    double perturber_inc = 80. * M_PI / 180.;
    reb_add_fmt(sim, "m a e inc", perturber_m, perturber_a, perturber_e, perturber_inc);

    // Add REBOUNDx effects
    // First tides_spin
    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);


    // Sun
    const double solar_k2 = 0.07;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.07 * star.m * star.r * star.r);

    const double solar_spin_period = 4.6 * 2. * M_PI / 365.;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin}); // Omega_x = Omega_y = 0 by default

    const double solar_Q = 1e6;
    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    // In the case of a spin that is synchronous with a circular orbit, tau is related to the tidal quality factor Q through the orbital mean motion n (see Lu et al. 2023 for discussion). Clearly that's not the case here, but gives us a reasonable starting point to start turning this knob
    double solar_tau = 1 / (2 * solar_Q * orb.n);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", solar_tau);

    // Planet
    const double planet_k2 = 0.4;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", planet_k2);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * planet_m * planet_r * planet_r);

    const double spin_period_p = 1. * 2. * M_PI / 365.; // days to reb years
    const double spin_p = (2. * M_PI) / spin_period_p;
    const double theta_p = 0. * M_PI / 180.;
    const double phi_p = 0. * M_PI / 180;
    struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_p, phi_p);
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_sv);

    const double planet_Q = 3e5;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1./(2.*planet_Q*orb.n));


    // add GR precession:
    struct rebx_force* gr = rebx_load_force(rebx, "gr_potential");
    rebx_add_force(rebx, gr);
    rebx_set_param_double(rebx, &gr->ap, "c", 10065.32); // in default units

    reb_move_to_com(sim);

    // Let's create a reb_rotation object that rotates to new axes with newz pointing along the total ang. momentum, and x along the line of
    // nodes with the invariable plane (along z cross newz)
    struct reb_vec3d newz = reb_vec3d_add(reb_tools_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)

    rebx_spin_initialize_ode(rebx, effect);

    system("rm -v output.txt");        // delete previous output file
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

      // orbits
      struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, *p1, *sun);
      double a1 = o1.a;
      double e1 = o1.e;
      double i1 = o1.inc;
      double Om1 = o1.Omega;
      double pom1 = o1.pomega;
      struct reb_vec3d n1 = o1.hvec;

      struct reb_particle com = reb_get_com_of_pair(sim->particles[0],sim->particles[1]);
      struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, *pert, com);
      double a2 = o2.a;
      double e2 = o2.e;
      double i2 = o2.inc;
      double Om2 = o2.Omega;
      double pom2 = o2.pomega;

      struct reb_vec3d* Omega_sun = rebx_get_param(rebx, sun->ap, "Omega");

      // Interpret planet spin in the rotating planet frame
      struct reb_vec3d* Omega_p_inv = rebx_get_param(rebx, p1->ap, "Omega");

      // Transform spin vector into planet frame, w/ z-axis aligned with orbit normal and x-axis aligned with line of nodes
      struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, n1);
      struct reb_rotation rot = reb_rotation_init_to_new_axes(n1, line_of_nodes); // Arguments to this function are the new z and x axes
      struct reb_vec3d srot = reb_vec3d_rotate(*Omega_p_inv, rot); // spin vector in the planet's frame

      // Interpret the spin axis in the more natural spherical coordinates
      double mag_p;
      double theta_p;
      double phi_p;
      reb_tools_xyz_to_spherical(srot, &mag_p, &theta_p, &phi_p);

      fprintf(of, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", sim->t, Omega_sun->x, Omega_sun->y, Omega_sun->z, mag_p, theta_p, phi_p, a1, e1, i1, Om1, pom1, a2, e2, i2, Om2, pom2); // print spins and orbits

      fclose(of);
    }

    if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
        reb_output_timing(sim, tmax);
    }
}
