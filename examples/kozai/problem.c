/**
 * Kozai cycles
 *
 * This example uses the IAS15 integrator to simulate
 * a Lidov Kozai cycle of a planet perturbed by a distant star.
 * The integrator automatically adjusts the timestep so that
 * even very high eccentricity encounters are resolved with high
 * accuracy.
 */
 #include <stdio.h>
 #include <stdlib.h>
 #include <unistd.h>
 #include <math.h>
 #include "rebound.h"
 #include "reboundx.h"
 #include "spin.c"

void heartbeat(struct reb_simulation* r);
double tmax = 1e5 * 2 * M_PI;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->dt             = M_PI*1e-2;     // initial timestep
    r->integrator        = REB_INTEGRATOR_IAS15;
    r->heartbeat        = heartbeat;

    // Initial conditions

    struct reb_particle star = {0};
    star.m  = 0.32;
    star.r = 0.5 * 0.00465;
    reb_add(r, star);

    // The planet (a zero mass test particle)
    double planet_m  = 0.05 * 9.55e-4; // in Jupiter masses
    double planet_r = 0.3 * 4.676e-4;
    double planet_a = 2.;
    double planet_e = 0.01;
    double planet_omega = 0. * (M_PI/180);
    reb_add_fmt(r, "m r a e omega", planet_m, planet_r, planet_a, planet_e, planet_omega);

    // The perturber
    struct reb_particle perturber = {0};
    double perturber_inc = 65. * (M_PI / 180.);
    double perturber_mass = 10. * 9.55e-4;
    double perturber_a  = 50.;
    double perturber_e = 0.52;
    reb_add_fmt(r, "m a e inc", perturber_mass, perturber_a, perturber_e, perturber_inc);

    struct rebx_extras* rebx = rebx_attach(r);

    struct rebx_force* effect = rebx_load_force(rebx, "spin");
    rebx_add_force(rebx, effect);
    // Sun
    const double solar_spin_period = 4.6 * 2. * M_PI / 365.;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    const double solar_k2 = 0.1;
    //const double solar_tau = 0.2 / solar_k2 * (2 * M_PI / 3.154e7); // seconds to reb years
    rebx_set_param_double(rebx, &r->particles[0].ap, "k2", solar_k2);
    // rebx_set_param_double(rebx, &r->particles[0].ap, "sigma", 6303.);
    rebx_set_param_double(rebx, &r->particles[0].ap, "moi", 0.07 * star.m * star.r * star.r);
    rebx_set_param_double(rebx, &r->particles[0].ap, "spin_sx", solar_spin * 0.0);
    rebx_set_param_double(rebx, &r->particles[0].ap, "spin_sy", solar_spin * 0.0);
    rebx_set_param_double(rebx, &r->particles[0].ap, "spin_sz", solar_spin * 1.0);
    // rebx_set_param_double(rebx, &r->particles[0].ap, "gr_spin", 1);
    //rebx_set_time_lag(r, rebx, &r->particles[0], solar_tau);
    rebx_set_q(r, rebx, &r->particles[0], &r->particles[1], 3e6);

    // P1
    const double spin_period_p = 1. * 2. * M_PI / 365.; // days to reb years
    const double spin_p = (2. * M_PI) / spin_period_p;
    const double planet_k2 = 0.4;
    //const double planet_tau = 0.02 / planet_k2 * (2 * M_PI / 3.154e7); // seconds to reb years
    const double planet_q = 3e5;
    const double theta_1 = 0. * M_PI / 180.; // initialize at one degree obliquity
    const double phi_1 = 0. * M_PI / 180;
    rebx_set_param_double(rebx, &r->particles[1].ap, "k2", planet_k2);
    //rebx_set_param_double(rebx, &r->particles[1].ap, "sigma", 1.75e15);
    rebx_set_param_double(rebx, &r->particles[1].ap, "moi", 0.25 * planet_m * planet_r * planet_r);
    rebx_set_param_double(rebx, &r->particles[1].ap, "spin_sx", spin_p * sin(theta_1) * sin(phi_1));
    rebx_set_param_double(rebx, &r->particles[1].ap, "spin_sy", spin_p * sin(theta_1) * cos(phi_1));
    rebx_set_param_double(rebx, &r->particles[1].ap, "spin_sz", spin_p * cos(theta_1));
    // rebx_set_time_lag(r, rebx, &r->particles[1], planet_tau);
    // rebx_set_param_double(rebx, &r->particles[1].ap, "gr_spin", 1);
    rebx_set_q(r, rebx, &r->particles[1], &r->particles[0], planet_q);


    // add GR:
    struct rebx_force* gr = rebx_load_force(rebx, "gr_full");
    rebx_add_force(rebx, gr);

    rebx_set_param_double(rebx, &gr->ap, "c", 10065.32); // in default units

    reb_move_to_com(r);
    rebx_align_simulation(r, rebx);
    rebx_spin_initialize_ode(r, effect);

    FILE* f = fopen("12_3_simple_nep_kozai_gr.txt","w");
    fprintf(f, "t,star_sx,star_sy,star_sz,magstar,a1,i1,e1,s1x,s1y,s1z,mag1,pom1,Om1,f1,a2,i2,e2,Om2,pom2\n");

    for (int i=0; i<5000000; i++){
        struct reb_particle* sun = &r->particles[0];
        struct reb_particle* p1 = &r->particles[1];
        struct reb_particle* pert = &r->particles[2];

        double* star_sx = rebx_get_param(rebx, sun->ap, "spin_sx");
        double* star_sy = rebx_get_param(rebx, sun->ap, "spin_sy");
        double* star_sz = rebx_get_param(rebx, sun->ap, "spin_sz");

        double* sx1 = rebx_get_param(rebx, p1->ap, "spin_sx");
        double* sy1 = rebx_get_param(rebx, p1->ap, "spin_sy");
        double* sz1 = rebx_get_param(rebx, p1->ap, "spin_sz");

        struct reb_orbit o1 = reb_tools_particle_to_orbit(r->G, *p1, *sun);
        double a1 = o1.a;
        double Om1 = o1.Omega;
        double i1 = o1.inc;
        double pom1 = o1.pomega;
        double f1 = o1.f;
        double e1 = o1.e;

        struct reb_vec3d s1 = {*sx1, *sy1, *sz1};

        struct reb_particle com = reb_get_com_of_pair(r->particles[0],r->particles[1]);
        struct reb_orbit o2 = reb_tools_particle_to_orbit(r->G, *pert, com);
        double a2 = o2.a;
        double Om2 = o2.Omega;
        double i2 = o2.inc;
        double pom2 = o2.pomega;
        double e2 = o2.e;

        // Interpret in the planet frame
        double magstar = sqrt((*star_sx) * (*star_sx) + (*star_sy) * (*star_sy) + (*star_sz) * (*star_sz));
        double mag1 = sqrt((*sx1) * (*sx1) + (*sy1) * (*sy1) + (*sz1) * (*sz1));
        double ob1 = acos(s1.z / mag1) * (180 / M_PI);

        if (i % 50000 == 0){
            printf("t=%e\t a1=%.6f\t o1=%0.5f\n", r->t / (2 * M_PI), a1, ob1);
            struct reb_vec3d gv = rebx_tools_spin_and_orbital_angular_momentum(r, rebx);
            double gv_mag = sqrt(gv.x*gv.x+gv.y*gv.y+gv.z*gv.z);
            printf("Tot orbital and spin ang mom: %0.10e %0.10e %0.10e\n", gv.x / gv_mag, gv.y / gv_mag, gv.z / gv_mag);
        }
        fprintf(f, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%.e,%e\n", r->t / (2 * M_PI), *star_sx, *star_sy, *star_sz, magstar, a1, i1, e1, s1.x, s1.y, s1.z, mag1, pom1, Om1, f1, a2, i2, e2, Om2, pom2);
        reb_integrate(r, r->t+(10 * 2 * M_PI));
    }
   rebx_free(rebx);
   reb_free_simulation(r);
}

void heartbeat(struct reb_simulation* r){
    if(reb_output_check(r, 20.*M_PI)){        // outputs to the screen
        // reb_output_timing(r, tmax);
    }
    /*
    if(reb_output_check(r, 12.)){            // outputs to a file
        reb_output_orbits(r, "orbits.txt");
    }
    */
}
