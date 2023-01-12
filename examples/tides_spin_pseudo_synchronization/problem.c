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
    const double p1_inc = 0.01;
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
    // The three components of the anguarl spin frequency (Omega_x, Omega_y, Omega_z)
    const double solar_k2 = 0.07;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    const double solar_spin_period = 27 * 2 * M_PI / 365; // 27 Days in REBOUND time units
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin}); // Omega.x = Omega.y = 0 by default

    // While not technically necessary, the fully dimensional moment of inertia (I) should also be set
    // This is required to evolve the spin axis! Without setting this value, the spin axis will remain stationary.
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.07 * solar_mass * solar_rad * solar_rad);

    // Finally, the last parameter is the constant time lag tau
    const double solar_Q = 1e6; // Often tidal dissipation is expressed in terms of a tidal quality factor Q
    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    // In the case of a spin that is synchronous with a circular orbit, tau is related to the tidal quality factor Q through the orbital mean motion n
    double solar_tau = 1 / (2 * solar_Q * orb.n);
    // See Lu et. al (2023) for discussion of the cases in which this approximation relating Q and tau is valid

    rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", solar_tau);

    // Planet - all of the above applies here too
    const double spin_period_1 = 0.5 * 2. * M_PI / 365.; // 0.5 days in reb years
    const double spin_1 = (2. * M_PI) / spin_period_1;
    const double planet_Q = 10000.;
    const double theta_1 = 30. * (M_PI / 180.);
    const double phi_1 = 0 * (M_PI / 180);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.3);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * p1_mass * p1_rad * p1_rad);

    // Let's consider a tilted planetary spin axis. The spin frequency vector is expressed in the same inertial reference frame as the simulation.
    // Given a polar angle theta (from the z axis) and azimuthal angle phi (from the x axis), we could either manually set spin axis components:
    struct reb_vec3d Omega_1;
    Omega_1.x = spin_1 * sin(theta_1) * cos(phi_1);
    Omega_1.y = spin_1 * sin(theta_1) * sin(phi_1);
    Omega_1.z = spin_1 * cos(theta_1);

    // Or use the built-in convenience function, which returns the Cartesian coordinates of the spin vector given:
    // magnitude, obliquity, and phase angle
    Omega_1 = reb_tools_spherical_to_xyz(spin_1, theta_1, phi_1);

    // For your own use case, you would just use whichever of the above two methods is most convenient (here we've done it three ways)
    // whichever of the two methods we use above, we need to remember to set the spin axis values
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_1);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1./(2*planet_Q*orb.n));

    reb_move_to_com(sim);

    // Let's create a reb_rotation object that rotates our simulation to new axes with newz pointing along the total ang. momentum, and x along the line of
    // nodes with the invariable plane (along z cross newz)
    struct reb_vec3d newz = reb_vec3d_add(reb_tools_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
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
      struct reb_vec3d* Omega_inv = rebx_get_param(rebx, p->ap, "Omega");

      // Transform spin vector into planet frame, w/ z-axis aligned with orbit normal and x-axis aligned with line of nodes
      struct reb_vec3d orbit_normal = orb.hvec;
      struct reb_vec3d line_of_nodes = reb_vec3d_cross((struct reb_vec3d){.z =1}, orbit_normal);
      struct reb_rotation rot = reb_rotation_init_to_new_axes(orbit_normal, line_of_nodes); // Arguments to this function are the new z and x axes
      struct reb_vec3d srot = reb_vec3d_rotate(*Omega_inv, rot); // spin vector in the planet's frame

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
