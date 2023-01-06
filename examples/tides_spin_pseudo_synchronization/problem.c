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
double tmax = 1000 * 2 * M_PI;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();

    // Star
    const double solar_mass = 1.;
    const double solar_rad = 0.00465; // Bodies with structure require radius! This is one solar radius
    reb_add_fmt(sim, "m r", solar_mass, solar_rad);// Central object

    // Fiducial hot Jupiter
    const double p1_mass = 1. * 9.55e-4; // in Jupiter masses * 1 Jupiter Mass / 1 Solar Mass
    const double p1_rad = 1. * 4.676e-4; // in Jupiter rad * 1 jupiter rad / 1 AU
    const double p1_e = 0.01;
    const double p1_inc = 0;
    reb_add_fmt(sim, "m a e inc r", p1_mass, 0.04072, p1_e, p1_inc, p1_rad); // Planet 1

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
    // The three components of spin (sx, sy, sz)
    const double solar_k2 = 0.07;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    const double solar_spin_period = 27 * 2 * M_PI / 365; // 27 Days in REBOUND time units
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sx", 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sy", 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sz", solar_spin);

    // While not technically necessary, the fully dimensional moment of inertia (moi) should also be set
    // This is required to evolve the spin axis! Without setting this value, the spin axis will remain stationary.
    rebx_set_param_double(rebx, &sim->particles[0].ap, "moi", 0.07 * solar_mass * solar_rad * solar_rad);

    // Finally, the last parameter is the tidal dissipation parameter (sigma)
    // This quantity is related to the constant time lag tau, and in the case of circular synchronous orbits, to the tital quality factor Q
    const double solar_Q = 1e6; // the tidal quality factor Q
    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    const double solar_tau = 1 / (2 * solar_Q * orb.n); // the constant time lag tau
    // See Lu et. al (2023) for discussion of the cases in which this approximation relating Q and tau is valid

    // You can calculate sigma directly (see Lu et al. (2023))
    double solar_sigma = 4 * sim->G * solar_tau / (3 * solar_rad * solar_rad * solar_rad * solar_rad * solar_rad * solar_k2);

    // Alternatively, you can use a helper function to calculate sigma from either Q or tau.
    // In both cases, you need to make sure that k2 and the physical radius are already set for the body. If they aren't you'll get a reb_error

    // To use rebx_tides_calc_sigma_from_Q, we need
    // 1) The rebx extras object
    // 2) The particle you wish to set sigma for
    // 3) The primary around which the body is orbiting
    // 4) The tidal quality factor Q you will use
    // Here we are recalculating our solar_sigma from above using the helper function (you would just choose most convenient option):
    solar_sigma = rebx_tides_calc_sigma_from_Q(rebx, &sim->particles[0], &sim->particles[1], solar_Q);

    // rebx_tides_calc_sigma_from_tau takes, as input:
    // 1) The rebx extras object
    // 2) The body you wish to set sigma for
    // 3) The Q value you will use
    solar_sigma = rebx_tides_calc_sigma_from_tau(rebx, &sim->particles[0], solar_tau);

    // For your own use case, you would just use whichever of the above three methods is most convenient (here we've done it three ways)
    // whichever of the three methods we use above, we need to remember to set sigma
    rebx_set_param_double(rebx, &sim->particles[0].ap, "sigma", solar_sigma);

    // Planet - all of the above applies here too
    const double spin_period_1 = 0.5 * 2. * M_PI / 365.; // 0.5 days in reb years
    const double spin_1 = (2. * M_PI) / spin_period_1;
    const double planet_Q = 10000.;
    const double theta_1 = 30. * (M_PI / 180.);
    const double phi_1 = 0 * (M_PI / 180);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.3);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "moi", 0.25 * p1_mass * p1_rad * p1_rad);

    // You can either manually set spin axis components:
    double planet_sx = spin_1 * sin(theta_1) * cos(phi_1);
    double planet_sy = spin_1 * sin(theta_1) * sin(phi_1);
    double planet_sz = spin_1 * cos(theta_1);

    // Or use the built-in convenience function, which returns the Cartesian coordinates of the spin vector given:
    // magnitude, obliquity, and phase angle
    struct reb_vec3d p_sv = reb_tools_spherical_to_xyz(spin_1, theta_1, phi_1);
    planet_sx = p_sv.x;
    planet_sy = p_sv.y;
    planet_sz = p_sv.z;

    // For your own use case, you would just use whichever of the above two methods is most convenient (here we've done it three ways)
    // whichever of the two methods we use above, we need to remember to set the spin axis values
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sx", planet_sx);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sy", planet_sy);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sz", planet_sz);

    double planet_sigma = rebx_tides_calc_sigma_from_Q(rebx, &sim->particles[1], &sim->particles[0], planet_Q);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "sigma", planet_sigma);

    reb_move_to_com(sim);

    // Let's create a reb_rotation object that rotates to new axes with newz pointing along the total ang. momentum, and x along the line of
    // nodes with the invariable plane (along z cross newz)
    struct reb_vec3d newz = rebx_tools_total_angular_momentum(rebx);
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)
    rebx_spin_initialize_ode(rebx, effect); // We must call this function before starting our integration to track the spins

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

      struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, *p, *star);
      double a = orb.a;
      double Om = orb.Omega;
      double inc = orb.inc;
      double pom = orb.pomega;
      double e = orb.e;

      // The spin vector in the inertial (in this case, invariant frame)
      double* sx_inv = rebx_get_param(rebx, p->ap, "sx");
      double* sy_inv = rebx_get_param(rebx, p->ap, "sy");
      double* sz_inv = rebx_get_param(rebx, p->ap, "sz");
      struct reb_vec3d spin_inv = {*sx_inv, *sy_inv, *sz_inv};

      // Transform spin vector into planet frame, w/ z-axis aligned with orbit normal and x-axis aligned with line of nodes
      struct reb_vec3d orbit_normal = orb.hvec;
      struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, orbit_normal);
      struct reb_rotation rot = reb_rotation_init_to_new_axes(orbit_normal, line_of_nodes); // Arguments to this function are the new z and x axes
      struct reb_vec3d srot = reb_vec3d_rotate(spin_inv, rot); // spin vector in the planet's frame

      // Interpret the spin axis in the more natural spherical coordinates
      double mag;
      double theta;
      double phi;
      reb_tools_xyz_to_spherical(srot, &mag, &theta, &phi);

      fprintf(of, "%e,%e,%e,%e,%e,%e,%e,%e,%e\n", sim->t, a, inc, e, pom, Om, mag, theta, phi);
      fclose(of);
    }

    if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
        reb_output_timing(sim, tmax);
    }
}
