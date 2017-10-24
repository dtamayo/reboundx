#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "rebound.h"
#include "reboundx.h"
#include "tools.h"
#define Rs 0.00464912633 // the solar radius in terms of AU
#define RJ 0.103*0.00464912633 // Jovian radius = 0.103Rs

double tmax = 5.0e9;

void heartbeat(struct reb_simulation* sim);
int merge(struct reb_simulation* const r, struct reb_collision c){
    printf("***hi***\n");
    return reb_collision_resolve_merge(r, c);
}

int main(int argc, char* argv[]){
   
    struct reb_simulation* sim = reb_create_simulation();
    srand(0); 
    sim->G = 4*M_PI*M_PI;
    sim->integrator = REB_INTEGRATOR_IAS15;
    //sim->ri_ias15.epsilon = 1e-10;
    sim->collision = REB_COLLISION_DIRECT;
    sim->collision_resolve = merge;//reb_collision_resolve_merge;
    sim->heartbeat = heartbeat;
    sim->dt = 1.e-3;
    
    
    reb_configure_box(sim,200,1,1,1);
    
    // put star on the center
    struct reb_particle star;
    star.m = 1.;//1.e8;
    star.r = Rs;
    reb_add(sim,star);
    struct reb_particle primary = sim->particles[0]; // using star as reference body(primary)
    
    // put n equal-mass planets(initial)
    
    double m = 10.*(1/332946.0487);// all initial planet masses 10. ME(fixed)
    double a[5],f[5]; // creating 4 SMA arrays and 4 phase arrays
    
    a[0] = 13.2,a[1]=32.3,a[2]=64.2,a[3]=73.7,a[4]=91.0; // initial semi-major axis in terms of AU
    int i;
    for(i=0;i < 5;i++){
        f[i] = reb_random_uniform(0.0,2.0*M_PI);// initial frequency in terms of redian for planets
    }
    
    double e = 0.0; // eccentricity for planet1
    double Omega = 0.0; // if orbits uniformly inclined, arbitary
    double omega = 0.0; // probably arbitary because initially circular orbit
    double inc= reb_random_uniform(0.0,M_PI/180.0); // inclination angle in terms of radian for planet1
    
    struct reb_particle p[5];
    for(i=0;i < 5;i++){
        p[i] = reb_tools_orbit_to_particle(sim->G, primary, m, a[i], e, inc, Omega, omega, f[i]);
        reb_add(sim,p[i]);
        p[i].r = RJ;
        p[i].hash = i+1;
        printf("%f\t%f \n",inc, f[i]);
        
    }
    
    reb_move_to_com(sim);
    /*for(int i=1;i<5;i++){
        struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[i], sim->particles[0]);
        printf("%d\t%f\t%f\t%f\n", i, o.a, o.e, o.inc);
    }*/
    /*struct rebx_extras* rebx = rebx_init(sim);
     
    double tau_d = 3.0e6; // 3.0Myr
    
    struct rebx_effect* paramsm = rebx_add(rebx,"modify_mass");
    double tau_m0 = 0.85e6; // 0.85Myr
     
    double* tau_m1= rebx_add_param(&sim->particles[1],"tau_mass",REBX_TYPE_DOUBLE);
    *tau_m1= tau_m0*(exp((sim->t)/tau_d));
    double* tau_m2= rebx_add_param(&sim->particles[2],"tau_mass",REBX_TYPE_DOUBLE);
    *tau_m2= tau_m0*(exp((sim->t)/tau_d));
    double* tau_m3= rebx_add_param(&sim->particles[3],"tau_mass",REBX_TYPE_DOUBLE);
    *tau_m3= tau_m0*(exp((sim->t)/tau_d));
    double* tau_m4= rebx_add_param(&sim->particles[4],"tau_mass",REBX_TYPE_DOUBLE);
    *tau_m4= tau_m0*(exp((sim->t)/tau_d));
    double* tau_m5= rebx_add_param(&sim->particles[5],"tau_mass",REBX_TYPE_DOUBLE);
    *tau_m5= tau_m0*(exp((sim->t)/tau_d));
    
    
    struct rebx_effect* paramse = rebx_add(rebx,"modify_orbits_forces");
    double tau_e0 = 1.0e6; // 1Myr
    
    double* tau_e1 = rebx_add_param(&sim->particles[1], "tau_e", REBX_TYPE_DOUBLE);
    *tau_e1= tau_e0*(exp((sim->t)/tau_d));
    double* tau_e2 = rebx_add_param(&sim->particles[2], "tau_e", REBX_TYPE_DOUBLE);
    *tau_e2= tau_e0*(exp((sim->t)/tau_d));
    double* tau_e3 = rebx_add_param(&sim->particles[3], "tau_e", REBX_TYPE_DOUBLE);
    *tau_e3= tau_e0*(exp((sim->t)/tau_d));
    double* tau_e4 = rebx_add_param(&sim->particles[4], "tau_e", REBX_TYPE_DOUBLE);
    *tau_e4= tau_e0*(exp((sim->t)/tau_d));
    double* tau_e5= rebx_add_param(&sim->particles[5], "tau_e", REBX_TYPE_DOUBLE);
    *tau_e5= tau_e0*(exp((sim->t)/tau_d));
    */ 
    //int*coordinates = rebx_add_param(paramse, "coordinates", REBX_TYPE_INT);
    //*coordinates = REBX_COORDINATES_PARTICLE;
    
    reb_integrate(sim, tmax);
    //rebx_free(rebx);
    reb_free_simulation(sim);
    return 0;
}


void heartbeat(struct reb_simulation* sim){
    if (reb_output_check(sim, 10000*sim->dt)){
        for(int i=1;i<5;i++){
            struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[i], sim->particles[0]);
            printf("%d\t%f\t%f\t%f\n", i, o.a, o.e, o.inc);
        }
        //reb_output_timing(sim, tmax);
    }
}
