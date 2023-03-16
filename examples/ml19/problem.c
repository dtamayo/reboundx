/**
 * Self-Consistent Spin-Tidal and Dynamical equations of motion (Eggleton et. al 1998)
 *
 * In particular, this simulates the obliquity sculpting of the Kepler multis due to convergent migration over 4 Myr.
 * Two planets are initialized just wide of the 3:2 MMR (0.173 and 0.233 AU) are migrated inward. After 2 Myr, migration is turned off. After a period of chaotic obliquity evolution, the inner planet is excited to high obliquity and maintained indefinitely
 * Result is based on Figure 3 in Millholland & Laughlin (2019), and reproduces Figure 3 in Lu et. al (2023)
 * For a more in-depth description of the various parameters that can be set in this simulation, please see the ipython examples for consistent tides & spin (any notebook with the prefix TidesSpin)
 * and migration (Migration.ipynb)
 *
 * integration time is artificially shortened to run quickly. Full run in paper with tmax = 4e6 orbits takes ~ 2 days
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "tides_spin.c"

void heartbeat(struct reb_simulation* sim);
//double tmax = 100 * 2 * M_PI; // set short to run quickly. Set to 4e6 * 2 * M_PI in paper
double tmax = 4e6 * 2 * M_PI;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Exact parameters from Millholland & Laughlin (2019)
    const double solar_mass = 1.;
    const double solar_rad = 0.00465;
    reb_add_fmt(sim, "m r", solar_mass, solar_rad);// Central object

    const double p1_mass = 5. * 3.0e-6; // in Earth masses * 1 Earth Mass / 1 Solar Mass
    const double p1_rad = 2.5 * 4.26e-5; // in Earth rad * 1 Earth rad / 1 AU
    reb_add_fmt(sim, "m a e r inc Omega pomega M", p1_mass, 0.17308688, 0.01, p1_rad, 0.5 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.)); // Planet 1

    const double p2_mass = 5. * 3.0e-6;
    const double p2_rad = 2.5 * 4.26e-5;
    reb_add_fmt(sim, "m a e r inc Omega pomega M", p2_mass, 0.23290608, 0.01, p2_rad, -0.431 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.)); // Planet 2
    sim->N_active = 3;
    sim->integrator = REB_INTEGRATOR_BS;
    sim->ri_bs.eps_abs = 1e-10;
    sim->ri_bs.eps_rel = 1e-10;
    sim->dt = 1e-3;
    sim->heartbeat = heartbeat;

    // Add REBOUNDx Additional effects
    // First Spin
    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);
    // Exact parameters from Millholland & Laughlin (2019)
    // Sun
    const double solar_spin_period = 20 * 2 * M_PI / 365;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    const double solar_Q = 1000000.;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", 0.07);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.07 * solar_mass * solar_rad * solar_rad);
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z = solar_spin}); // Omega_x = Omega_y = 0 by default

    // We assume tau = 1/(2*n*Q) with n the mean motion, even though the spin is not synchronized with the orbit (see Lu et al. (2023))
    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", 1./(2.*orb.n*solar_Q));

    // P1
    const double spin_period_1 = 5. * 2. * M_PI / 365.; // 5 days in REBOUND time units
    const double spin_1 = (2. * M_PI) / spin_period_1;
    const double planet_Q = 10000.;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * p1_mass * p1_rad * p1_rad);

    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1./(2.*orb.n*planet_Q));
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", (struct reb_vec3d){.y=spin_1 * -0.0261769, .z=spin_1 * 0.99965732});

    // P2
    double spin_period_2 = 3. * 2. * M_PI / 365.; // 3 days in REBOUND time units
    double spin_2 = (2. * M_PI) / spin_period_2;
    rebx_set_param_double(rebx, &sim->particles[2].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "I", 0.25 * p2_mass * p2_rad * p2_rad);
    rebx_set_param_vec3d(rebx, &sim->particles[2].ap, "Omega", (struct reb_vec3d){.y=spin_2 * 0.0249736, .z=spin_2 * 0.99968811});

    struct reb_orbit orb2 = reb_tools_particle_to_orbit(sim->G, sim->particles[2], sim->particles[0]);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau", 1./(2.*orb2.n*planet_Q));
/*
    double mag2 = spin_2;//7.3e1;
    double theta2 = 1.745411e-02;
    double phi2 = 0.0;//-1.570796;

    // Define the inverse transformation: planet frame -> inv frame right now
    struct reb_vec3d lon2 = reb_vec3d_cross((struct reb_vec3d){.z =1}, orb2.hvec);  // Line of nodes is the new x-axis
    struct reb_rotation rot2 = reb_rotation_init_to_new_axes(orb2.hvec, lon2);      // Arguments to this function are the new z and x axes
    struct reb_rotation rot_test = reb_rotation_inverse(rot2);

    struct reb_vec3d temp_p = reb_tools_spherical_to_xyz(mag2, theta2, phi2);
    struct reb_vec3d sv2 = reb_vec3d_rotate(temp_p, rot_test);
    rebx_set_param_vec3d(rebx, &sim->particles[2].ap, "Omega", sv2);
    */

    // And migration
    struct rebx_force* mo = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mo);

    // Set migration parameters
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -5e6 * 2 * M_PI);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", (-5e6 * 2 * M_PI) / 1.1);

    reb_move_to_com(sim);

    // Let's create a reb_rotation object that rotates to new axes with newz pointing along the total ang. momentum, and x along the line of
    // nodes with the invariable plane (along z cross newz)
    struct reb_vec3d newz = reb_vec3d_add(reb_tools_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    rebx_simulation_irotate(rebx, rot); // This rotates our simulation into the invariable plane aligned with the total ang. momentum (including spin)
    rebx_spin_initialize_ode(rebx, effect);

    // Run simulation
    //system("rm -v output_orbits.txt"); // remove previous output files
    //system("rm -v output_spins.txt");

    // Simulation Archive
    reb_simulationarchive_automate_interval(sim,"archive_01_25_low_eps.bin",1000.);

    reb_integrate(sim, tmax/2);

    printf("Migration Switching Off\n");
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", INFINITY);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", INFINITY);

    reb_integrate(sim, tmax);

    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
  if(reb_output_check(sim, tmax/100000)){        // outputs every 100 REBOUND years
    struct rebx_extras* const rebx = sim->extras;
    FILE* of_orb = fopen("output_orbits_01_25_BS_low_eps.txt", "a");
    FILE* of_spins = fopen("output_spins_01_25_BS_low_eps.txt", "a");
    if (of_orb == NULL || of_spins == NULL){
        reb_error(sim, "Can not open file.");
        return;
    }

    struct reb_particle* sun = &sim->particles[0];
    struct reb_particle* p1 = &sim->particles[1];
    struct reb_particle* p2 = &sim->particles[2];

    // Orbit information
    struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, *p1, *sun);
    double a1 = o1.a;
    double e1 = o1.e;
    double i1 = o1.inc;
    double Om1 = o1.Omega;
    double pom1 = o1.pomega;
    struct reb_vec3d norm1 = o1.hvec;

    struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, *p2, *sun);
    double a2 = o2.a;
    double e2 = o2.e;
    double i2 = o2.inc;
    double Om2 = o2.Omega;
    double pom2 = o2.pomega;
    struct reb_vec3d norm2 = o2.hvec;

    // Spin vectors - all initially in invariant plane
    struct reb_vec3d* Omega_sun = rebx_get_param(rebx, sun->ap, "Omega");

    // Interpret both planet spin vectors in the rotating planet frame in spherical coordinates
    struct reb_vec3d* Omega_p1 = rebx_get_param(rebx, p1->ap, "Omega");

    struct reb_vec3d lon1 = reb_vec3d_cross((struct reb_vec3d){.z =1}, norm1);  // Line of nodes is the new x-axis
    struct reb_rotation rot1 = reb_rotation_init_to_new_axes(norm1, lon1);      // Arguments to this function are the new z and x axes
    struct reb_vec3d sv1 = reb_vec3d_rotate(*Omega_p1, rot1);

    double mag1;
    double theta1;
    double phi1;
    reb_tools_xyz_to_spherical(sv1, &mag1, &theta1, &phi1);

    struct reb_vec3d* Omega_p2 = rebx_get_param(rebx, p2->ap, "Omega");

    struct reb_vec3d lon2 = reb_vec3d_cross((struct reb_vec3d){.z =1}, norm2); // Line of nodes is the new x-axis
    struct reb_rotation rot2 = reb_rotation_init_to_new_axes(norm2, lon2); // Arguments to this function are the new z and x axes
    struct reb_vec3d sv2 = reb_vec3d_rotate(*Omega_p2, rot2);

    double mag2;
    double theta2;
    double phi2;
    reb_tools_xyz_to_spherical(sv2, &mag2, &theta2, &phi2);

    fprintf(of_orb, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", sim->t, a1, e1, i1, pom1, Om1, norm1.x, norm1.y, norm1.z, a2, e2, i2, pom2, Om2, norm2.x, norm2.y, norm2.z);  // prints the spins and orbits of all bodies
    fclose(of_orb);
    fprintf(of_spins, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", sim->t, Omega_sun->x, Omega_sun->y, Omega_sun->z, mag1, theta1, phi1, mag2, theta2, phi2, Omega_p1->x, Omega_p1->y, Omega_p1->z, Omega_p2->x, Omega_p2->y, Omega_p2->z);
    fclose(of_spins);
  }

  //if(reb_output_check(sim, 100.*M_PI)){        // outputs to the screen
  //    reb_output_timing(sim, tmax);
  //}
}
