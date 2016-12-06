/**
 * Restarting simulations
 * 
 * This example demonstrates how to restart a simulation
 */
#include <stdio.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "core.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    struct reb_particle p = {0}; 
    p.m     = 1.;   
    reb_add(sim, p); 

    double m = 0.;
    double a1 = 1.;
    double a2 = 2.;
    double e = 0.4;
    double omega = 0.;
    double f = 0.;

    struct reb_particle p1 = reb_tools_orbit2d_to_particle(sim->G, p, m, a1, e, omega, f);
    struct reb_particle p2 = reb_tools_orbit2d_to_particle(sim->G, p, m, a2, e, omega, f);
    reb_add(sim,p1);
    reb_add(sim,p2);
    reb_move_to_com(sim);

    struct rebx_extras* rebx = rebx_init(sim);
    
    struct rebx_effect* effect1 = rebx_add(rebx, "modify_orbits_direct");    					// directly update particles' orbital elements each timestep
	
    double* tau_a = rebx_add_param(effect1, "tau_a", REBX_TYPE_DOUBLE); 			// add semimajor axis damping on inner planet (e-folding timescale)
	double* tau_omega = rebx_add_param(effect1, "tau_omega", REBX_TYPE_DOUBLE);	// add linear precession (set precession period). Won't do anything for modify_orbits_forces
	double* tau_e = rebx_add_param(effect1, "tau_e", REBX_TYPE_DOUBLE);			// add eccentricity damping on particles[2] (e-folding timescale)

    *tau_a = 3.;
    *tau_omega = 4.;
    *tau_e = 5.;
    
    struct rebx_effect* effect2 = rebx_add(rebx, "modify_orbits_forces");
    double* a = rebx_add_param(effect2, "a", REBX_TYPE_DOUBLE);
    double* b = rebx_add_param(effect2, "b", REBX_TYPE_DOUBLE);
    double* c = rebx_add_param(effect2, "c", REBX_TYPE_DOUBLE);
    double* d = rebx_add_param(effect2, "d", REBX_TYPE_DOUBLE);
    
    *a = 100.;
    *b = 200.;
    *c = 300.;
    *d = 400.;
    
    struct rebx_effect* effect3 = rebx_add(rebx, "gr");
    int length = 3;
    int* integ = rebx_add_param_(effect3, "integ", REBX_TYPE_INT, 1, &length);
    
    int shape2D[2] = {2,3};
    double* doub = rebx_add_param_(effect3, "doub", REBX_TYPE_DOUBLE, 2, shape2D);
    
    int shape3D[3] = {2,4,3};
    uint32_t* uinteg = rebx_add_param_(effect3, "uinteg", REBX_TYPE_UINT32, 3, shape3D);
    
    struct reb_orbit* orbit = rebx_add_param_(effect3, "orbit", REBX_TYPE_ORBIT, 1, &length);
    
    int integs[3] = {-10,0,10};
    memcpy(integ, integs, sizeof(*integ)*3);
    
    double doubs[6] = {1.,2.,3.,4.,5.,6.};
    memcpy(doub, doubs, sizeof(*doub)*6);
    
    uint32_t uintegs[24] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
    memcpy(uinteg, uintegs, sizeof(*uinteg)*24);
    
    
    struct reb_orbit orbits[3] = {{.a=1., .e=0.1}, {.a=2., .e=0.2}, {.a=3., .e=0.3}};
    memcpy(orbit, orbits, sizeof(*orbit)*3);
    
    struct reb_particle* planet = &sim->particles[2];
    
    double* ap = rebx_add_param(planet, "ap", REBX_TYPE_DOUBLE);
    *ap = 100.;
    int* integ2 = rebx_add_param_(planet, "integ2", REBX_TYPE_INT, 1, &length);
    double* doub2 = rebx_add_param_(planet, "doub2", REBX_TYPE_DOUBLE, 2, shape2D);
    uint32_t* uinteg2 = rebx_add_param_(planet, "uinteg2", REBX_TYPE_UINT32, 3, shape3D);
    struct reb_orbit* orbit2 = rebx_add_param_(planet, "orbit2", REBX_TYPE_ORBIT, 1, &length);
    memcpy(integ2, integs, sizeof(*integ)*3);
    memcpy(doub2, doubs, sizeof(*doub)*6);
    memcpy(uinteg2, uintegs, sizeof(*uinteg)*24);
    memcpy(orbit2, orbits, sizeof(*orbit)*3);
    
    //reb_integrate(sim,20.);
    rebx_output_binary(rebx, "test.bin");
    
    rebx_free(rebx);
    reb_free_simulation(sim);
    
    sim = reb_create_simulation();
    p.m     = 1.;
    reb_add(sim,p);
    reb_add(sim,p1);
    reb_add(sim,p2);
    reb_move_to_com(sim);
    planet = &sim->particles[2];
    
    rebx = rebx_create_extras_from_binary(sim, "test.bin");
    
    effect1 = rebx_get_effect(rebx, "modify_orbits_direct");
    tau_a = rebx_get_param(effect1, "tau_a");
    tau_omega = rebx_get_param(effect1, "tau_omega");
    tau_e = rebx_get_param(effect1, "tau_e");
    
    effect3 = rebx_get_effect(rebx, "gr");
    integ = rebx_get_param(effect3, "integ");
    doub = rebx_get_param(effect3, "doub");
    uinteg = rebx_get_param(effect3, "uinteg");
    orbit = rebx_get_param(effect3, "orbit");
    
    integ2 = rebx_get_param(planet, "integ2");
    doub2 = rebx_get_param(planet, "doub2");
    uinteg2 = rebx_get_param(planet, "uinteg2");
    orbit2 = rebx_get_param(planet, "orbit2");
    ap = rebx_get_param(planet, "ap");
    
    printf("%f\t%f\t%f\n", *tau_a, *tau_omega, *tau_e);
    printf("\n");

    for(int i=0; i<3; i++){
        printf("%d\t", integ[i]);
    }
    printf("\n");
    
    for(int i=0; i<6; i++){
        printf("%f\t", doub[i]);
    }
    printf("\n");
    
    for(int i=0; i<24; i++){
        printf("%u\t", uinteg[i]);
    }
    printf("\n");
    
    for(int i=0; i<3; i++){
        printf("%f\t", orbit[i].e);
    }
    printf("\n");
    
    printf("*******\n");
    
    printf("%f\n", *ap);
    for(int i=0; i<3; i++){
        printf("%d\t", integ2[i]);
    }
    printf("\n");
    
    for(int i=0; i<6; i++){
        printf("%f\t", doub2[i]);
    }
    printf("\n");
    
    for(int i=0; i<24; i++){
        printf("%u\t", uinteg2[i]);
    }
    printf("\n");
    
    for(int i=0; i<3; i++){
        printf("%f\t", orbit2[i].e);
    }
    printf("\n");
    
    struct rebx_param* param = rebx_get_param_node(planet, "integ2");
    
    printf("%s\n", param->name);
    printf("%d\t%d\n", param->shape[0], param->shape[1]);
    
    struct rebx_param* param2 = rebx_get_param_node(planet, "doub2");
    printf("%d\t%d\n", param2->shape[0], param2->shape[1]);
}
