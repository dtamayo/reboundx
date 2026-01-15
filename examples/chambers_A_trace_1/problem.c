#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"


char TITLE[100] = "chambers_A_";

double get_radii(double m, double rho){
    return pow((3*m)/(4*M_PI*rho),1./3.);
}

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();

    r->integrator = REB_INTEGRATOR_TRACE;
    r->ri_trace.S_peri = reb_integrator_trace_switch_peri_none;
    r->ri_trace.r_crit_hill = 3.*1.21;
    r->collision = REB_COLLISION_DIRECT;
    r->G = 39.476926421373;
    r->dt = 6./365.;
    r->rand_seed = 1;

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
    double min_frag_mass = 1.4e-8;
    rebx_set_param_double(rebx, &fragmenting->ap, "fc_min_frag_mass", min_frag_mass);
    rebx_set_param_pointer(rebx, &fragmenting->ap, "fc_particle_list_file", "family_tree.csv");


    //Assigning mass and number of planetary embryos and planetesimals
    double rho = 5.05e6; //3 g/cm^3
    int n_emb = 14; //number of planetary embryos
    double m_emb = 2.8e-7; //mass of each embryo
    int n_pl = 140; //number of planetesimals
    double m_pl = 2.8e-08;
    int ef = 1; //expansion factor, in case we want to use any
    double a_pl[n_pl];  // Array to hold the values of planetesimal semi-major axis
    double a_emb[n_emb]; // Same for embryos

    FILE *file_pl = fopen("original_chambers_2013_pl_semis.txt", "r");
    for (int i = 0; i < n_pl; i++) {
        fscanf(file_pl, "%lf", &a_pl[i]);}
    // Close the file after reading
    fclose(file_pl);

    FILE *file_emb = fopen("original_chambers_2013_emb_semis.txt", "r");
    for (int i = 0; i < n_emb; i++) {
        fscanf(file_emb, "%lf", &a_emb[i]);}
    // Close the file after reading
    fclose(file_emb);


    struct reb_particle star = {0};
    star.m = 1.00;
    star.r = 0.00465;
    reb_simulation_add(r, star);
    
    // Add planetary embryos
    
    for (int i=0; i<n_emb; i++){
        double a = a_emb[i];      // semi major axis
        double m = m_emb;
        double inc = reb_random_uniform(r,0,0.0175);
        double ecc = reb_random_uniform(r,0,0.01);
        double omega = reb_random_uniform(r,0,2*M_PI);
        double Omega = reb_random_uniform(r,0,2*M_PI);
        double f = reb_random_uniform(r,0,2*M_PI);
        //now build particle from orbit
        struct reb_particle emb = reb_particle_from_orbit(r->G, star, m, a, ecc, inc, Omega, omega, f);
        emb.r = get_radii(m, rho)*ef;
        reb_simulation_add(r, emb); 
    }

    //add planetesimals
    for (int i=0; i<n_pl; i++){
        double a = a_pl[i];      // semi major axis
        double m = m_pl;
        double inc = reb_random_uniform(r,0,0.0175);
        double ecc = reb_random_uniform(r,0,0.01);
        double omega = reb_random_uniform(r,0,2*M_PI);
        double Omega = reb_random_uniform(r,0,2*M_PI);
        double f = reb_random_uniform(r,0,2*M_PI);
        //now build particle from orbit
        struct reb_particle pl = reb_particle_from_orbit(r->G, star, m, a, ecc, inc, Omega, omega, f);
        pl.r = get_radii(m, rho)*ef;
        
        reb_simulation_add(r, pl);
    }

    
    double m,a,e,inc,Omega,omega,f; //Omega=longitude of ascending node, omega= argument of pericenter in RADIANS

    //Add Jupiter and Saturn
    struct reb_particle Jup = reb_particle_from_orbit(r->G, r->particles[0], m=9.543e-4, a=5.20349, e=0.048381, inc=0.365*(M_PI/180), Omega=0.0, omega=68.3155*(M_PI/180), f=227.0537*(M_PI/180));
    Jup.r = get_radii(Jup.m, rho); 
    Jup.hash = reb_hash("JUPITER");
    reb_simulation_add(r,Jup);

    struct reb_particle Sat = reb_particle_from_orbit(r->G, r->particles[0], m=0.0002857, a=9.54309, e=0.052519, inc=0.8892*(M_PI/180), Omega=M_PI, omega=324.5263*(M_PI/180),f=256.9188*(M_PI/180));
    Sat.r = get_radii(Sat.m, rho);
    Sat.hash = reb_hash("SATURN");
    reb_simulation_add(r,Sat);


    reb_simulation_move_to_com(r);  // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    double run_time = 10;
    reb_simulation_save_to_file_interval(r,TITLE,1.e5);
    reb_simulation_integrate(r, run_time);
}