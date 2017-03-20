/**
 * Saving and loading simulations
 * 
 * This example demonstrates how to restart a simulation with all REBOUNDx effects and parameters.
 */
#include <stdio.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "core.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation(); // make a simple sim with star and 1 planet
    struct reb_particle p = {0}; 
    p.m     = 1.;   
    reb_add(sim, p);
    struct reb_particle p1 = {0};
    p1.x=1.;
    reb_add(sim, p1);
    struct reb_particle* planet = &sim->particles[1];
    struct rebx_extras* rebx = rebx_init(sim);
   
    // Add some (arbitrary) parameters to effects and particles
    struct rebx_effect* gr = rebx_add(rebx, "gr");
    double* a = rebx_add_param(gr, "a", REBX_TYPE_DOUBLE);
    double* b = rebx_add_param(gr, "b", REBX_TYPE_DOUBLE);
    int* c = rebx_add_param(gr, "c", REBX_TYPE_INT);
    *a = 1.;
    *b = 2.;
    *c = 3;
    
    printf("a: Original = %f\n", *a);
    printf("b: Original = %f\n", *b);
    printf("c: Original = %d\n", *c);

    // See array_parameters example for details on adding array parameters
    int ndim = 3;
    int shape3D[3] = {2,4,3};
    int* integers = rebx_add_param_array(planet, "integers", REBX_TYPE_INT, ndim, shape3D);
    int ints[24] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
    memcpy(integers, ints, sizeof(*integers)*24);

    rebx->integrator = REBX_INTEGRATOR_RK4;
    gr->force_as_operator = 1;
    gr->operator_order = 2;

    printf("rebx->integrator: Original = %d\n", rebx->integrator);
    printf("gr->force_as_operator: Original = %d\n", gr->force_as_operator);
    printf("gr->operator_order: Original = %d\n", gr->operator_order);
    
    // We now have to save both a REBOUND binary (for the simulation) and a REBOUNDx one (for parameters and effects)
    reb_output_binary(sim, "reb.bin");
    rebx_output_binary(rebx, "rebx.bin");
    
    rebx_free(rebx);
    reb_free_simulation(sim);
   
    // We now reload the simulation and the rebx instance (which adds previously loaded effects to the simulation)
    sim = reb_create_simulation_from_binary("reb.bin");
    rebx = rebx_create_extras_from_binary(sim, "rebx.bin");
    
    gr = rebx_get_effect(rebx, "gr");
    a = rebx_get_param(gr, "a");
    b = rebx_get_param(gr, "b");
    c = rebx_get_param(gr, "c");
    
    printf("rebx->integrator: Loaded = %d\n", rebx->integrator);
    printf("gr->force_as_operator: Loaded = %d\n", gr->force_as_operator);
    printf("gr->operator_order: Loaded = %d\n", gr->operator_order);
    
    planet = &sim->particles[1];
    integers = rebx_get_param(planet, "integers");
    
    printf("a: Loaded = %f\n", *a);
    printf("b: Loaded = %f\n", *b);
    printf("c: Loaded = %d\n", *c);
    
    printf("Original planet integers:\t");
    for(int i=0; i<24; i++){
        printf("%d\t", ints[i]);
    }
    printf("\nLoaded planet integers:\t\t");
    for(int i=0; i<24; i++){
        printf("%d\t", integers[i]);
    }
    printf("\n");

    // You would now integrate as usual
}
