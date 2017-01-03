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
    struct reb_simulation* sim = reb_create_simulation(); // we begin by adding 3 arbitrary particles to a simulation
    struct reb_particle p = {0}; 
    p.m     = 1.;   
    reb_add(sim, p);
    struct reb_particle p1 = {0};
    p1.x=1.;
    reb_add(sim, p1);
    struct reb_particle p2 = {0};
    p2.x=1.;
    reb_add(sim, p2);
    
    struct rebx_extras* rebx = rebx_init(sim);
    
    struct rebx_effect* gr = rebx_add(rebx, "gr");
    double* a = rebx_add_param(gr, "a", REBX_TYPE_DOUBLE);
    double* b = rebx_add_param(gr, "b", REBX_TYPE_DOUBLE);
    double* c = rebx_add_param(gr, "c", REBX_TYPE_DOUBLE);
    *a = 1.;
    *b = 2.;
    *c = 3.;
    
    printf("a: Original = %f\n", *a);
    printf("b: Original = %f\n", *b);
    printf("c: Original = %f\n", *c);
    
    struct reb_particle* planet = &sim->particles[2];
    int shape1D = 3;
    int ints[3] = {-10,0,10};
    int* integers = rebx_add_param_(planet, "integers", REBX_TYPE_INT, 1, &shape1D);
    memcpy(integers, ints, sizeof(*integers)*3);
    
    struct reb_orbit* orbits = rebx_add_param_(planet, "orbits", REBX_TYPE_ORBIT, 1, &shape1D);
    struct reb_orbit orbs[3] = {{.a=1., .e=0.1}, {.a=2., .e=0.2}, {.a=3., .e=0.3}};
    memcpy(orbits, orbs, sizeof(*orbits)*3);
    
    int shape2D[2] = {2,3};
    double doubs[6] = {1.,2.,3.,4.,5.,6.};
    double* doubles = rebx_add_param_(planet, "doubles", REBX_TYPE_DOUBLE, 2, shape2D);
    memcpy(doubles, doubs, sizeof(*doubles)*6);
    
    int shape3D[3] = {2,4,3};
    uint32_t* uintegers = rebx_add_param_(planet, "uintegers", REBX_TYPE_UINT32, 3, shape3D);
    uint32_t uints[24] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
    memcpy(uintegers, uints, sizeof(*uintegers)*24);
  
    reb_output_binary(sim, "reb.bin");
    rebx_output_binary(rebx, "rebx.bin");
    
    rebx_free(rebx);
    reb_free_simulation(sim);
    
    sim = reb_create_simulation_from_binary("reb.bin");
    rebx = rebx_create_extras_from_binary(sim, "rebx.bin");
    
    gr = rebx_get_effect(rebx, "gr");
    a = rebx_get_param(gr, "a");
    b = rebx_get_param(gr, "b");
    c = rebx_get_param(gr, "c");
    
    planet = &sim->particles[2];
    integers = rebx_get_param(planet, "integers");
    doubles = rebx_get_param(planet, "doubles");
    uintegers = rebx_get_param(planet, "uintegers");
    orbits = rebx_get_param(planet, "orbits");
    
    printf("a: Loaded = %f\n", *a);
    printf("b: Loaded = %f\n", *b);
    printf("c: Loaded = %f\n", *c);
    
    printf("Originals:\t");
    for(int i=0; i<3; i++){
        printf("%d\t", ints[i]);
    }
    printf("\nLoaded:\t\t");
    for(int i=0; i<3; i++){
        printf("%d\t", integers[i]);
    }
    printf("\n");
    
    printf("Originals:\t");
    for(int i=0; i<6; i++){
        printf("%f\t", doubs[i]);
    }
    printf("\nLoaded:\t\t");
    for(int i=0; i<6; i++){
        printf("%f\t", doubles[i]);
    }
    printf("\n");
    
    printf("Originals:\t");
    for(int i=0; i<24; i++){
        printf("%d\t", uints[i]);
    }
    printf("\nLoaded:\t\t");
    for(int i=0; i<24; i++){
        printf("%d\t", uintegers[i]);
    }
    printf("\n");
}
