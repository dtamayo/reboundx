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


void compute_transformation_angles(struct reb_simulation* sim, double* theta1, double* theta2){
    // From celmech line 330
    struct reb_vec3d gtot_vec = reb_tools_angular_momentum(sim);
    double gtot = sqrt(gtot_vec.x * gtot_vec.x + gtot_vec.y * gtot_vec.y + gtot_vec.z * gtot_vec.z);
    double ghat_x = gtot_vec.x / gtot;
    double ghat_y = gtot_vec.y / gtot;
    double ghat_z = gtot_vec.z / gtot;
    double ghat_perp = sqrt(1 - ghat_z * ghat_z);
    *theta1 = M_PI / 2 - atan2(ghat_y, ghat_x);
    *theta2 = M_PI / 2 - atan2(ghat_z, ghat_perp);
}

struct reb_vec3d EulerAnglesTransform(struct reb_vec3d xyz, double Omega, double I, double omega){
    // celmech line 341
    double x = xyz.x;
    double y = xyz.y;
    double z = xyz.z;

    double s1 = sin(omega);
    double c1 = cos(omega);
    double x1 = c1 * x - s1 * y;
    double y1 = s1 * x + c1 * y;
    double z1 = z;

    double s2 = sin(I);
    double c2 = cos(I);
    double x2 = x1;
    double y2 = c2 * y1 - s2 * z1;
    double z2 = s2 * y1 + c2 * z1;

    double s3 = sin(Omega);
    double c3 = cos(Omega);
    double x3 = c3 * x2 - s3 * y2;
    double y3 = s3 * x2 + c3 * y2;
    double z3 = z2;

    struct reb_vec3d shifted = {x3, y3, z3};
    return shifted;
}

void align_simulation(struct reb_simulation* sim){
    // celmech line 360
    const int N_real = sim->N - sim->N_var;
    double theta1, theta2;
    compute_transformation_angles(sim, &theta1, &theta2);
    for (int i = 0; i <= N_real; i++){
        struct reb_particle* p = &sim->particles[i];
	struct reb_vec3d pos = {p->x, p->y, p->z};
	struct reb_vec3d vel = {p->vx, p->vy, p->vz};
  }
}

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
    const double solar_mass = 1.;
    const double solar_rad = 0.00465;
    reb_add_fmt(sim, "m r", solar_mass, solar_rad);// Central object

    const double p1_mass = 5. * 3.0e-6; // in Earth masses * 1 Earth Mass / 1 Solar Mass
    const double p1_rad = 2.5 * 4.26e-5; // in Earth rad * 1 Earth rad / 1 AU
    reb_add_fmt(sim, "m a e r inc Omega pomega M", p1_mass, 0.17308688, 0.01, p1_rad, 0.5 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.)); // Planet 1

    const double p2_mass = 5. * 3.0e-6;
    const double p2_rad = 2.5 * 4.26e-5;
    reb_add_fmt(sim, "m a e r inc Omega pomega M", p2_mass, 0.23290608, 0.01, p2_rad, -0.431 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.)); // Planet 2
    reb_move_to_com(sim);
    sim->N_active = 3;
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->dt = 1e-2;

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
    rebx_set_param_double(rebx, &sim->particles[0].ap, "moi", 0.07 * solar_mass * solar_rad * solar_rad);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "spin_sx", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "spin_sy", solar_spin * 0.0);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "spin_sz", solar_spin * 1.0);

    // P1
    double spin_period_1 = 5. * 2. * M_PI / 365.; // 5 days in reb years
    double spin_1 = (2. * M_PI) / spin_period_1;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "q", 10000.);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "moi", 0.25 * p1_mass * p1_rad * p1_rad);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "spin_sx", spin_1 * 0.0);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "spin_sy", spin_1 * -0.0261769);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "spin_sz", spin_1 * 0.99965732);

    // P2
    double spin_period_2 = 3. * 2. * M_PI / 365.; // 3 days in reb years
    double spin_2 = (2. * M_PI) / spin_period_2;
    rebx_set_param_double(rebx, &sim->particles[2].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "q", 10000.);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "moi", 0.25 * p2_mass * p2_rad * p2_rad);
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

    FILE* f = fopen("8_11_sm_test.txt","w");
    fprintf(f, "t,a1,i1,e1,sx1,sy1,sz1,S1,pomega1,Omega1,f1,x1,y1,z1,a2,i2,e2,sx2,sy2,sz2,S2,pomega2,Omega2,f2,x2,y2,z2\n");
    int cond = 0;
     for (int i=0; i<100000; i++){

         struct reb_particle* sun = &sim->particles[0];
         struct reb_particle* p1 = &sim->particles[1];
         struct reb_particle* p2 = &sim->particles[2];

	 if (sim->t > 2e6 * 2 * M_PI && cond == 0){
		printf("Migration Switching Off\n");
		rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", INFINITY);
    		rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", INFINITY);
	 	cond = 1;
	 }

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

	 //printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", sim->t, a1, i1, e1, s1_inv.x, s1_inv.y, s1_inv.z, sqrt(s1_inv.x * s1_inv.x + s1_inv.y * s1_inv.y + s1_inv.z * s1_inv.z), pom1, Om1, p1->x, p1->y, p1->z, a2, i2, e2, s2_inv.x, s2_inv.y, s2_inv.z, sqrt(s2_inv.x * s2_inv.x + s2_inv.y * s2_inv.y + s2_inv.z * s2_inv.z), pom2, Om2, p2->x, p2->y, p2->z);


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
    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* r){
    if(reb_output_check(r, 25)){        // outputs to the screen
        //reb_output_timing(r, 1e4);
    }
}
