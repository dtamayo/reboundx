/**
 * Array Parameters 
 * 
 * This example shows how to create, populate and retrieve array parameters.
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
    struct reb_particle p1 = {0};
    p1.x=1.;
    reb_add(sim, p1);
    struct reb_particle* planet = &sim->particles[1];
   
    struct rebx_extras* rebx = rebx_init(sim);
    
    // To add arrays we have to additionally pass the number of dimensions in the array and its shape. First 1D:  
    int ndim = 1;
    int shape1D = 3;
    int ints[3] = {-10,0,10};
    // In addition to specifying the object to attach the parameter to, the name, and type, we pass ndim and shape
    int* integers = rebx_add_param_array(planet, "integers", REBX_TYPE_INT, ndim, &shape1D);
    // REBOUNDx allocates and manages the memory for all parameters, so need to copy the array to the returned memory location
    memcpy(integers, ints, sizeof(*integers)*3);
    
    // integers is a pointer to the parameter's contents, so changing it will also update the parameter 
    integers[0] = -1000;

    // If we ever need to get the pointer again, we call get_param
    integers = rebx_get_param(planet, "integers");

    printf("integers:\t");
    for(int i=0; i<3; i++){
        printf("%d ", integers[i]);
    }
    printf("\n");

    // We can also add arrays of any variable type defined in the rebx_param_type enum (reboundx.h)
    struct reb_orbit* orbits = rebx_add_param_array(planet, "orbits", REBX_TYPE_ORBIT, 1, &shape1D);
    struct reb_orbit orbs[3] = {{.a=1., .e=0.1}, {.a=2., .e=0.2}, {.a=3., .e=0.3}};
    memcpy(orbits, orbs, sizeof(*orbits)*3);
   
    // For multidimensional arrays, we have to pass a shape array with the dimensions (in row major order)
    // Multidimensional arrays are stored as unwrapped 1D arrays, so have to be entered this way.
    // Here we add a 2 row by 3 col array of [[1,2,3], [4,5,6]]
    int shape2D[2] = {2,3};
    double doubs[6] = {1.,2.,3.,4.,5.,6.};
    double* doubles = rebx_add_param_array(planet, "doubles", REBX_TYPE_DOUBLE, 2, shape2D);
    memcpy(doubles, doubs, sizeof(*doubles)*6);

    // When you add a parameter, REBOUNDx automatically creates metadata for it, which can be useful for arrays. 
    // You can get the parameter with all its metadata with
    struct rebx_param* doub_param = rebx_get_param_node(planet, "doubles");

    // We can now use the array's strides to access the array more intuitively
    int* strides = doub_param->strides;
    // e.g. we can get the [1,2] element (6, since all indices start at 0)
    printf("doubles element [1,2] = %f\n", doubles[1*strides[0] + 2*strides[1]]);
   
    // This becomes more useful with higher dimensions
    int shape3D[3] = {2,4,3};
    uint32_t* uintegers = rebx_add_param_array(planet, "uintegers", REBX_TYPE_UINT32, 3, shape3D);
    uint32_t uints[24] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
    memcpy(uintegers, uints, sizeof(*uintegers)*24);
    
    struct rebx_param* uints_param = rebx_get_param_node(planet, "uintegers");
    strides = uints_param->strides;
    printf("uintegers element [1,2,1] = %u\n", uints[1*strides[0] + 2*strides[1] + 1*strides[2]]);
}
