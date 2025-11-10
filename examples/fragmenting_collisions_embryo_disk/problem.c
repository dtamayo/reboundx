#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"


char TITLE[100] = "embryo_disk_";

double get_radii(double m, double rho){
    return pow((3*m)/(4*M_PI*rho),1./3.);
}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 10.)){
        //reb_simulation_output_timing(r, 0);
        printf("Walltime(s) = %f %d\n", r->walltime, r->N);
    }
}

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    r->integrator = REB_INTEGRATOR_MERCURIUS;
    r->collision = REB_COLLISION_DIRECT;
    r->G = 39.476926421373;
    r->dt = 6./365.;

    // The random seed is passed as a command line argument
    if (argc == 2){
        r->rand_seed = atoi(argv[1]);
        strcat(TITLE, argv[1]);
    }

    //Setting the collision resolve module
    struct rebx_extras* rebx = rebx_attach(r);
    struct rebx_collision_resolve* fragmenting = rebx_load_collision_resolve(rebx, "fragmenting_collisions");
    rebx_add_collision_resolve(rebx, fragmenting);

    //Choose minimum fragment mass
    rebx_set_param_double(rebx, &fragmenting->ap, "fc_min_frag_mass", 0.01);

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

    double rho = 5.05e6; //3 g/cm^3

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
        emb.r = get_radii(m, rho);

        reb_simulation_add(r, emb);
    }

    //Optional: initiate the family tree file, to save initial particle IDs as well
    FILE* of = fopen("family_tree.csv", "a");
    fprintf(of, "particle_id, parent_t, parent_p,\n"); 
    for(int i=0; i<r->N; i++){
        struct reb_particle* p = &(r->particles[i]); // First object in collision
        rebx_fragmenting_collisions_set_new_id(r, fragmenting, &r->particles[i]);
        int* new_id = rebx_get_param(rebx, p->ap, "fc_id");
        fprintf(of, "%d,", *new_id);
        fprintf(of, " , ");
        fprintf(of, " , ");
        fprintf(of, "\n"); 
    }
    fclose(of);

    reb_simulation_move_to_com(r);  // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    double run_time = 1e5;
    reb_simulation_save_to_file_interval(r,TITLE,1.e2);
    reb_simulation_integrate(r, run_time);


    reb_simulation_free(r);
}
