#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"


char TITLE[100] = "embryo_trace_";
char TITLE_FAMTREE[100] = "famtree_";
char TITLE_HEARTBEAT[100] = "heartbeat_";
double RHO = 5.05e6; //3 g/cm^3

double get_radii(double m, double rho){
    return pow((3*m)/(4*M_PI*rho),1./3.);
}

void heartbeat(struct reb_simulation* sim){
    if (reb_simulation_output_check(sim, 10)){
        FILE *fp = fopen("heartbeat_output.txt", "a");  // append mode
        if (fp == NULL){
            perror("Error opening file");
            return;
        }
        struct reb_vec3d angular_momentum = reb_simulation_angular_momentum(sim);
        double total_mass = 0;
        for (int i = 0; i < sim->N; i++){
            total_mass = total_mass + sim->particles[i].m;
        }
        fprintf(fp, "%e,", sim->t);
        fprintf(fp, "%f,", sim->walltime);
        fprintf(fp, "%d,", sim->N);
        fprintf(fp, "%e,", total_mass);
        fprintf(fp, "%e, ", angular_momentum.x);
        fprintf(fp, "%e, ", angular_momentum.y);
        fprintf(fp, "%e, ", angular_momentum.z);
        fprintf(fp, "%e\n", reb_simulation_energy(sim));

        fclose(fp);
    }
    
}

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    // The random seed is passed as a command line argument
    if (argc == 2){
        r->rand_seed = atoi(argv[1]);
        strcat(TITLE, argv[1]);
	strcat(TITLE_FAMTREE, argv[1]);
    }

    r->integrator = REB_INTEGRATOR_TRACE;
    r->ri_trace.S_peri = reb_integrator_trace_switch_peri_none;
    r->ri_trace.r_crit_hill = 3.*1.21;
    r->collision = REB_COLLISION_DIRECT;
    r->G = 39.476926421373;
    r->dt = 6./365.;
    r->exact_finish_time = 0;
    //r->heartbeat = heartbeat;

    //Setting the collision resolve module
    struct rebx_extras* rebx = rebx_attach(r);
    struct rebx_collision_resolve* fragmenting = rebx_load_collision_resolve(rebx, "fragmenting_collisions");
    rebx_add_collision_resolve(rebx, fragmenting);

    //Choose minimum fragment mass
    double min_frag_mass = 1.4e-8;
    rebx_set_param_double(rebx, &fragmenting->ap, "fc_min_frag_mass", min_frag_mass);
    rebx_set_param_pointer(rebx, &fragmenting->ap, "fc_particle_list_file", TITLE_FAMTREE);


    //Assigning mass and number of planetary embryos and planetesimals
    struct reb_particle star = {0};
    star.m = 1.0;
    star.r = 0.00465;
    reb_simulation_add(r, star);

    // Constants for mass range
    double lunar_mass = 3.8e-8;
    double earth_mass = 3e-6;
    double mass_min = 0.6 * lunar_mass;
    double mass_max = 0.2 * earth_mass;

    // Add 30 planetary embryos
    for (int i = 0; i < 30; i++) {
        double a = reb_random_uniform(r, 0.1, 0.5);     // semi-major axis in AU
        double e = reb_random_uniform(r, 0.0, 0.01);    // eccentricity
        double inc = reb_random_uniform(r, 0.0, M_PI/180.);                        // inclination
        double omega = reb_random_uniform(r, 0.0, 2.*M_PI); // argument of periapsis
        double Omega = reb_random_uniform(r, 0.0, 2.*M_PI); // longitude of ascending node
        double f = reb_random_uniform(r, 0.0, 2.*M_PI);     // mean anomaly
        double m = reb_random_uniform(r, mass_min, mass_max);  // in solar masses

        struct reb_particle emb = reb_particle_from_orbit(r->G, star, m, a, e, inc, Omega, omega, f);
        emb.r = get_radii(m, RHO);

        reb_simulation_add(r, emb);
    }


    reb_simulation_move_to_com(r);  // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    double run_time = 10e5;
    reb_simulation_save_to_file_step(r,TITLE,50);
    reb_simulation_integrate(r, run_time);

    reb_simulation_free(r);
}
