/**
 * Constant time lag model for tides (Hut 1981)
 *
 * In particular, this simulates post-main sequence tidal interactions between the Earth and Sun near its tip-RGB phase.
 * Definitely see the corresponding ipython example, as well as the documentation, for more explanations along the way of the various parameters and assumptions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* sim);

struct reb_vec3d transform(double inc, double omega, struct reb_vec3d spin_inv){
    // This ts a vector from the INVARIANT frame to the PLANET frame
    double sx = spin_inv.x;
    double sy = spin_inv.y;
    double sz = spin_inv.z;

    double t[3][3];

    t[0][0] = cos(omega);
    t[0][1] = sin(omega);
    t[0][2] = 0;
    t[1][0] = -cos(inc) * sin(omega);
    t[1][1] = cos(inc) * cos(omega);
    t[1][2] = sin(inc);
    t[2][0] = sin(inc) * sin(omega);
    t[2][1] = -sin(inc) * cos(omega);
    t[2][2] = cos(inc);

    struct reb_vec3d spin_planet = {0};

    spin_planet.x = sx * t[0][0] + sy * t[0][1] + sz * t[0][2];
    spin_planet.y = sx * t[1][0] + sy * t[1][1] + sz * t[1][2];
    spin_planet.z = sx * t[2][0] + sy * t[2][1] + sz * t[2][2];

    return spin_planet;
}

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();

    // Units of AU, 2pi and Msun: G = 1
    // sim->G = 4*M_PI*M_PI;           // Units of AU, yr and Msun
    //sim->heartbeat = heartbeat;

    struct reb_particle sun = {0};
    sun.m = 1.0;
    sun.r = 0.00465;
    reb_add(sim, sun);

    struct reb_orbit p1o = {0};      // p1
    p1o.a = 0.173086;                     // in AU
    p1o.e = 0.01;
    p1o.inc = 0.5 * (M_PI / 180.);
    p1o.Omega = 0.0;
    p1o.pomega = 0.0;
    p1o.omega = p1o.Omega + p1o.pomega;
    p1o.M = 0;
    p1o.f = 0; //reb_tools_M_to_f(p1o.e, p1o.M);
    double p1_mass = 5. * 3.0e-6;
    struct reb_particle p1 = reb_tools_orbit_to_particle(sim->G, sun, p1_mass, p1o.a, p1o.e, p1o.inc, p1o.Omega, p1o.omega, p1o.f);
    p1.r = 2.5 * 4.26e-5;
    reb_add(sim, p1);

    struct reb_orbit p2o = {0};      // p1
    p2o.a = 0.23290608;                     // in AU
    p2o.e = 0.01;
    p2o.inc = 0.431 * (M_PI / 180.);
    p2o.Omega = 0.0;
    p2o.pomega = 0.0;
    p2o.omega = p2o.Omega + p2o.pomega;
    p2o.M = 0;
    p2o.f = 0; // reb_tools_M_to_f(p2o.e, p2o.M);
    double p2_mass = 5. * 3.0e-6;
    struct reb_particle p2 = reb_tools_orbit_to_particle(sim->G, sun, p2_mass, p2o.a, p2o.e, p2o.inc, p2o.Omega, p2o.omega, p2o.f);
    p2.r = 2.5 * 4.26e-5;
    reb_add(sim, p2);
    reb_move_to_com(sim);
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->dt = 1e-2;
    sim->N_active = 3;

    // Add REBOUNDx Additional effects
    // First Spin
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* effect = rebx_load_force(rebx, "spin");
    rebx_add_force(rebx, effect);
    // Sun
    double solar_spin_period = 20 * 2 * M_PI / 365;
    double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", 0.07);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "q", 100000.);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "moi", 0.07 * sun.m * sun.r * sun.r);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "spin_sx", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "spin_sy", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "spin_sz", solar_spin * 1.0);

    // P1
    double spin_period_1 = 5 * 2 * M_PI / 365; // 5 days in reb years
    double spin_1 = (2 * M_PI) / spin_period_1;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "q", 10000.);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "moi", 0.25 * p1_mass * p1.r * p1.r);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "spin_sx", spin_1 * 0.0);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "spin_sy", spin_1 * -0.0261769);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "spin_sz", spin_1 * 0.99965732);

    // P2
    double spin_period_2 = 3 * 2 * M_PI / 365; // 5 days in reb years
    double spin_2 = (2 * M_PI) / spin_period_2;
    rebx_set_param_double(rebx, &sim->particles[2].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "q", 10000.);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "moi", 0.25 * p2_mass * p2.r * p2.r);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "spin_sx", spin_2 * 0.0);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "spin_sy", spin_2 * 0.0249736);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "spin_sz", spin_2 * 0.99968811);

    // And migration
    struct rebx_force* mo = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mo);

    // Set migration parameters
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -5e6 * 2 * M_PI);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", (-5e6 * 2 * M_PI) / 1.1);

    // Run simulation
    rebx_spin_initialize_ode(sim, effect);
    struct reb_ode* spin_ode = rebx_get_param(rebx, effect->ap, "ode");
    printf("%f\n", spin_ode->length);
/*
    FILE* f = fopen("7_6_rebx_test_exact_sm19.txt","w");
    fprintf(f, "t,a1,i1,e1,sx1,sy1,sz1,S1,pomega1,Omega1,f1,x1,y1,z1,a2,i2,e2,sx2,sy2,sz2,S2,pomega2,Omega2,f2,x2,y2,z2\n");
     for (int i=0; i<100000; i++){

         struct reb_particle* sun = &sim->particles[0];
         struct reb_particle* p1 = &sim->particles[1];
         struct reb_particle* p2 = &sim->particles[2];

         double* sx1 = rebx_get_param(rebx, p1->ap, "spin_sx");
         double* sy1 = rebx_get_param(rebx, p1->ap, "spin_sy");
         double* sz1 = rebx_get_param(rebx, p1->ap, "spin_sz");

         double* sx2 = rebx_get_param(rebx, p2->ap, "spin_sx");
         double* sy2 = rebx_get_param(rebx, p2->ap, "spin_sy");
         double* sz2 = rebx_get_param(rebx, p2->ap, "spin_sz");

         struct reb_orbit o1 = reb_tools_particle_to_orbit(sim->G, *p1, *sun);
         double a1 = o1.a;//vis_viva(r, &p1, &sun);
         double Om1 = o1.Omega;
         double i1 = o1.inc;
         double pom1 = o1.pomega;
         double f1 = o1.f;
         double e1 = o1.e;

         struct reb_orbit o2 = reb_tools_particle_to_orbit(sim->G, *p2, *sun);
         double a2 = o2.a;//vis_viva(r, &p2, &sun);
         double Om2 = o2.Omega;
         double i2 = o2.inc;
         double pom2 = o2.pomega;
         double f2 = o2.f;
         double e2 = o2.e;

         struct reb_vec3d s1_inv = {*sx1, *sy1, *sz1};
         struct reb_vec3d s2_inv = {*sx2, *sy2, *sz2};

         struct reb_vec3d s1 = transform(i1, Om1, s1_inv);
         struct reb_vec3d s2 = transform(i2, Om2, s2_inv);

         // Interpret in the planet frame
         double mag1 = sqrt(s1.x * s1.x + s1.y * s1.y + s1.z * s1.z);
         double ob1 = acos(s1.z / mag1) * (180 / M_PI);
         double mag2 = sqrt(s2.x * s2.x + s2.y * s2.y + s2.z * s2.z);
         double ob2 = acos(s2.z / mag2) * (180 / M_PI);

         if (i % 100 == 0){
             printf("t=%f\t a1=%.6f\t a2=%.6f\t o1=%0.5f\t o2=%0.5f\n", sim->t / (2 * M_PI), a1, a2, ob1, ob2);
         }
         fprintf(f, "%f,%f,%f,%f,%0.10f,%0.10f,%0.10f,%0.10f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%0.10f,%0.10f,%0.10f,%0.10f,%f,%f,%f,%f,%f,%f\n", sim->t / (2 * M_PI), a1, i1, e1, s1.x, s1.y, s1.z, mag1, pom1, Om1, f1, p1->x,p1->y,p1->z,a2, i2, e2, s2.x, s2.y, s2.z, mag2, pom2, Om2, f2, p2->x,p2->y,p2->z);
         reb_integrate(sim,sim->t+(40 * 2 * M_PI));
     }
     */
    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* r){
    if(reb_output_check(r, 25)){        // outputs to the screen
        //reb_output_timing(r, 1e4);
    }
}
