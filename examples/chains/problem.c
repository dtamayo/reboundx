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
double tmax = 1e5 * 2 * M_PI;
double delta = 1.05;

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();

    // Star: M-dwarf
    const double solar_mass = 0.0898;
    const double solar_rad = 0.117 * 0.00465; // Bodies with structure require radius! This is one solar radius
    reb_add_fmt(sim, "m r", solar_mass, solar_rad);// Central object

    const double innermost_P = 9.206690 * 2. * M_PI / 365.; // f
    const double innermost_a = 3.705241e-2;
    //const double innermost_P = 12.35294 * 2. * M_PI / 365;

    // Trappist-1b
    const double b_mass = 1.374 * 3e-6;
    const double b_rad = 1.116 * 4.264e-5;
    const double b_a = 1.154e-2 * 1.1;
    const double b_e = 0.007;
    const double b_inc = (90. - 89.728) * (M_PI / 180.);
    reb_add_fmt(sim, "m r a e inc", b_mass, b_rad, b_a, b_e, b_inc);

    // Trappist-1c
    const double c_mass = 1.308 * 3e-6;
    const double c_rad = 1.097 * 4.264e-5;
    const double c_a = b_a * (8./5.) * delta;
    // const double c_a = 1.580e-2;
    const double c_e = 0.002;
    const double c_inc = (90. - 89.778) * (M_PI / 180.);
    reb_add_fmt(sim, "m r a e inc", c_mass, c_rad, c_a, c_e, c_inc);

    // Trappist-1d
    const double d_mass = 0.388 * 3e-6;
    const double d_rad = 0.788 * 4.264e-5;
    const double d_a = c_a * (5./3.) * delta;
    // const double d_a = 2.227e-2;
    const double d_e = 0.007;
    const double d_inc = (90. - 89.896) * (M_PI / 180.);
    reb_add_fmt(sim, "m r a e inc", d_mass, d_rad, d_a, d_e, d_inc);

    // Trappist-1e
    const double e_mass = 0.692 * 3e-6;
    const double e_rad = 0.920 * 4.264e-5;
    const double e_a = d_a * (3./2.) * delta;
    // const double e_a = 2.925e-2;
    const double e_e = 0.005;
    const double e_inc = (90. - 89.793) * (M_PI / 180.);
    reb_add_fmt(sim, "m r a e inc", e_mass, e_rad, e_a, e_e, e_inc);

    // Trappist-1f
    const double f_mass = 1.039 * 3e-6;
    const double f_rad = 1.045 * 4.264e-5;
    // const double f_P = innermost_P;
    const double f_a = e_a * (3./2.) * delta;
    // const double f_a = 3.849e-2;
    const double f_e = 0.01;
    const double f_inc = 2.0 * (M_PI / 180.);
    // const double f_inc = (90. - 89.740) * (M_PI / 180.);
    reb_add_fmt(sim, "m r a e inc", f_mass, f_rad, f_a, f_e, f_inc); // Planet 1


    // Trappist-1g
    const double g_mass = 1.321 * 3e-6;
    const double g_rad = 1.129 * 4.264e-5;
    // const double g_P = innermost_P * (3./2.) * delta;
    const double g_a = f_a * (4./3.) * delta;
    // const double g_a = 4.683e-2;//f_a * cbrt(((3. / 2.) * delta) * ((3. / 2.) * delta));
    const double g_e = 0.003;
    const double g_inc = -1.0 * (M_PI / 180.);//
    // const double g_inc = (90. - 89.742) * (M_PI / 180.);
    reb_add_fmt(sim, "m r a e inc", g_mass, g_rad, g_a, g_e, g_inc); // Planet 1

    // Trappist-1h
    const double h_mass = 0.326 * 3e-6;
    const double h_rad = 0.755 * 4.264e-5;
    // const double h_P = g_P * (3./2.) * delta;
    const double h_a = g_a * (3./2.) * delta;
    // const double h_a = 6.189e-2;//g_a * cbrt(((3. / 2.) * delta) * ((3. / 2.) * delta));
    const double h_e = 0.004;
    const double h_inc = (90. - 89.805) * (M_PI / 180.);//0.195 * (M_PI / 180.);
    reb_add_fmt(sim, "m r a e inc", h_mass, h_rad, h_a, h_e, h_inc);

    struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    double ts = orb.P / 15.12345; // timestep as a function of orbital period of innermost planet

    //double inner_edge = innermost_a - 0.0005;
    //double ts = sqrt(4. * M_PI / (solar_mass) * (inner_edge * inner_edge * inner_edge)) / 15.12345; // timestep as a function of period at inner disk edge

    sim->N_active = 8;
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->dt = ts;
    sim->heartbeat = heartbeat;

    // Add tides_spin as a REBOUNDx additional effect
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    // Star
    const double solar_k2 = 0.307;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);

    const double solar_spin_period = 3.3 * 2 * M_PI / 365;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin});

    const double rg_star = 0.2;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", rg_star * solar_mass * solar_rad * solar_rad);
    //rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", solar_tau);

    // Planets - all assumed to have same tidal parameters
    const double planet_k2 = 0.299;
    const double planet_tau = (712.37 / 86400) * 2 * M_PI / 365; // seconds->days->reb years
    const double planet_rg = 0.3308;

    for (int i = 1; i < sim->N_active; i++){
      rebx_set_param_double(rebx, &sim->particles[i].ap, "k2", planet_k2);
      rebx_set_param_double(rebx, &sim->particles[i].ap, "tau", planet_tau);
      rebx_set_param_double(rebx, &sim->particles[i].ap, "I", planet_rg * sim->particles[i].m * sim->particles[i].r * sim->particles[i].r);

      // tidally locked & 0 obliquity
      struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, sim->particles[i], sim->particles[0]);
      const double spin_period = orb.P;
      const double spin_rate = (2. * M_PI / spin_period);
      // printf("%e\n", orb.a);

      struct reb_vec3d norm = orb.hvec;
      struct reb_vec3d nhat = reb_vec3d_normalize(norm);
      // struct reb_vec3d Omega_p = reb_vec3d_mul(nhat, spin_rate);

      double breakup = sqrt(sim->G * sim->particles[i].m / (sim->particles[i].r * sim->particles[i].r * sim->particles[i].r));
      struct reb_vec3d Omega_p = reb_vec3d_mul(nhat, breakup);
      rebx_set_param_vec3d(rebx, &sim->particles[i].ap, "Omega", Omega_p);
    }


    // Migration - MODIFY_ORBITS FORCES
    // *****************************************************************************************************

    struct rebx_force* mo = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mo);

    double mig_rate = -3e4 * 2 * M_PI;
    double step = 0.1;
    double K = 125.;
    /*
    for (int i = 1; i < sim->N_active; i++){
      rebx_set_param_double(rebx, &sim->particles[i].ap, "tau_a", mig_rate);
      rebx_set_param_double(rebx, &sim->particles[i].ap, "tau_e", mig_rate / K);
      mig_rate /= (1 + i * step);
    }
    */

    // migration applied only to outer planet
    rebx_set_param_double(rebx, &sim->particles[7].ap, "tau_a", mig_rate);
    rebx_set_param_double(rebx, &sim->particles[7].ap, "tau_e", mig_rate / K);
    // *****************************************************************************************************


    // Migration - TYPE_I_MIGRATION
    // *****************************************************************************************************
    /*
    struct rebx_force* mig = rebx_load_force(rebx, "type_I_migration");

    // Using fiducial values for now
    //double inner_edge = innermost_a - 0.0005;

    rebx_set_param_double(rebx, &mig->ap, "ide_position", inner_edge);
    rebx_set_param_double(rebx, &mig->ap, "ide_width", 0.005);
    rebx_set_param_double(rebx, &mig->ap, "tIm_flaring_index", 0.0);
    rebx_set_param_double(rebx, &mig->ap, "tIm_surface_density_exponent", 1.5);
    rebx_set_param_double(rebx, &mig->ap, "tIm_surface_density_1", 0.00011255);
    rebx_set_param_double(rebx, &mig->ap, "tIm_scale_height_1", 0.01);
    rebx_add_force(rebx, mig);

    // double ts = sqrt(4. * M_PI / (solar_mass) * (inner_edge * inner_edge * inner_edge)) / 15.12345; // timestep as a function of period at inner disk edge
    */
    // *****************************************************************************************************

    reb_move_to_com(sim);

    // Rotate into invariant plane
    struct reb_vec3d newz = reb_vec3d_add(reb_tools_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);

    rebx_simulation_irotate(rebx, rot);
    rebx_spin_initialize_ode(rebx, effect);

    system("rm -v test.txt"); // remove previous output file
    FILE* of = fopen("test.txt", "a");
    //fprintf(of, "t,a1,i1,pom1,Om1,mag1,theta1,phi1,hx1,hx2,hx3,ex1,ey1,ez1,a2,i2,pom2,Om2,mag2,theta2,phi2,hx2,hy2,hz2,ex2,ey2,ez2,a3,i3,pom3,Om3,mag3,theta3,phi3,hx3,hy3,hz3,ex3,ey3,ez3\n");
    fprintf(of, "t");
    //for (int i = 1; i < sim->N_active; i++){
    //  fprintf(of, ",a%d,i%d,pom%d,Om%d,mag%d,theta%d,phi%d,hx%d,hy%d,hz%d,ex%d,ey%d,ez%d", i, i, i, i, i, i, i, i, i, i, i, i, i);
    //}
    //fprintf(of, "\n");
    fprintf(of, ",hx,hy,hz,sx,sy,sz\n");
    fclose(of);
    /*
    reb_integrate(sim, tmax/2);

    printf("\nMigration Switching Off\n");
    for (int i = 1; i < sim->N_active; i++){
      rebx_set_param_double(rebx, &sim->particles[i].ap, "tau_a", INFINITY);
      rebx_set_param_double(rebx, &sim->particles[i].ap, "tau_e", INFINITY);
    }
    */
    //reb_integrate(sim, tmax/10);

    // turn on migration
    //rebx_set_param_double(rebx, &sim->particles[3].ap, "tau_a", mig_rate);
    //rebx_set_param_double(rebx, &sim->particles[3].ap, "tau_e", mig_rate / K);

    reb_integrate(sim, tmax);
    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    // Output spin and orbital information to file
    if(reb_output_check(sim, 1)){        // outputs every REBOUND year
      struct rebx_extras* const rebx = sim->extras;
      FILE* of = fopen("test.txt", "a");
      if (of==NULL){
          reb_error(sim, "Can not open file.");
          return;
      }

      struct reb_particle* star = &sim->particles[0];
      fprintf(of, "%e", sim->t);
      for (int i = 1; i < sim->N_active; i++){
        struct reb_particle* p = &sim->particles[i];

        struct reb_orbit orb = reb_tools_particle_to_orbit(sim->G, *p, *star);
        double a = orb.a;
        double Om = orb.Omega;
        double inc = orb.inc;
        double om = orb.omega;
        struct reb_vec3d hvec = orb.hvec;
        struct reb_vec3d evec = orb.evec;

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
        //fprintf(of, ",%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e", a, inc, om, Om, mag, theta, phi, hvec.x, hvec.y, hvec.z, evec.x, evec.y, evec.z);
        //fprintf(of, ",%e,%e,%e,%e,%e", a, inc, e, pom, Om);
      }

      struct reb_vec3d angmom = reb_tools_angular_momentum(sim);
      struct reb_vec3d spinmom = rebx_tools_spin_angular_momentum(rebx);
      fprintf(of, ",%e,%e,%e,%e,%e,%e\n", angmom.x, angmom.y, angmom.x, spinmom.x, spinmom.y, spinmom.z);
      fclose(of);
    }

    if(reb_output_check(sim, 20.*M_PI)){        // outputs to the screen
        reb_output_timing(sim, tmax);
    }
}
