/**
 * Adding parameters to objects in REBOUNDx
 * 
 * This example walks through adding parameters to particles, forces and operators in REBOUNDx
 */
#include <stdio.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]){
    // We start by making a simulation and adding GR with REBOUNDx
    struct reb_simulation* sim = reb_create_simulation();

    struct reb_particle star = {0};
    star.m     = 1.;   
    reb_add(sim, star);

    struct reb_particle planet = {0};
    planet.x = 1.;
    planet.vy = 1.;

    reb_add(sim, planet);
    
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* gr = rebx_load_force(rebx, "gr");
    rebx_add_force(rebx, gr);
    
    /* The documentation page https://reboundx.readthedocs.io/en/latest/effects.html lists the various required and optional parameters that need to be set for each effect in REBOUNDx. 
     *
     * Parameters are stored in particles, forces and operators through the linked list 'ap', which is common to all three objects. 
     *
     * For simple types (int and double) we add them through their corresponding setters.
     * In both cases, we pass the rebx struct, a pointer to the head of the linked list that we want to modify, and the value we want to set.
     */

    double c = 10064.915; // speed of light in default units of AU, Msun and yr/2pi
    rebx_set_param_double(rebx, &gr->ap, "c", c); 

    // After setting the parameters we want to set, we would integrate as usual. 

    double tmax = 10.;
    reb_integrate(sim, tmax); 

    /* At any point, we can access the parameters we set (e.g., some effects could update these values as the simulation progresses). We get all parameter types back with rebx_get_param, which returns a void pointer that we are responsible for casting to the correct type. Since we are not modifying the linked list, we don't pass a reference to ap like above*/

    double* new_c = rebx_get_param(rebx, gr->ap, "c");

    /* It is important to check returned parameters!
     * If the parameter is not found, it will return a NULL pointer.
     * If this happens and you dereference it, you will get a segmentation fault.
     * If you get a seg fault, the first thing you should check are returned pointers.
     */

    if (new_c != NULL){
        printf("c=%f\n", *new_c);
    }

    /* The above functionality is probably enough for most users.
     * If you find you need to do more complicated things, read below!
     *
     * More complicated types can be added with set_param_pointer.
     * Instead of passing a value, we now pass a pointer to the variable.
     * As a simple example adding the gr pointer to particles[1] (with no meaning whatsoever):
     */

    rebx_set_param_pointer(rebx, &sim->particles[1].ap, "force", gr);

    /* In order to be able to write binary files, and for interoperability with Python, REBOUNDx keeps a registered list of parameter names and their corresponding types.
     * Above we have used parameter names that have been registered by various effects.
     * If you try to use a name that's not on the list, REBOUNDx will print an error and continue execution without adding the parameter.
     */

    rebx_set_param_int(rebx, &gr->ap, "q", 7);

    int* q = rebx_get_param(rebx, gr->ap, "q");

    if (q != NULL){
        printf("q=%d\n", *q);
    }
    else{
        printf("q is NULL because we didn't register the name first.\n");
    }

    /* You could register the name permanently in rebx_register_params in core.c, or you can do it manually in problem.c, passing a name and rebx_param_type enum (defined in reboundx.h):
     */

    rebx_register_param(rebx, "q", REBX_TYPE_INT);

    rebx_set_param_int(rebx, &gr->ap, "q", 7);
    q = rebx_get_param(rebx, gr->ap, "q");

    if (q != NULL){
        printf("q=%d\n", *q);
    }
    else{
        printf("q is NULL because we didn't register the name first.\n");
    }

    /* Finally, you might want to add your own custom structs (e.g. from another library).
     * This can be done straightforwardly by adding it as a pointer.
     */
    
    struct SPH_sim{
        double dt;
        int Nparticles;
    };

    struct SPH_sim my_sph_sim = {.dt = 0.1, .Nparticles=10000};

    rebx_register_param(rebx, "sph", REBX_TYPE_POINTER);
    rebx_set_param_pointer(rebx, &gr->ap, "sph", &my_sph_sim);

    // If you need to get it back:
    
    struct SPH_sim* sph = rebx_get_param(rebx, gr->ap, "sph");

    if (sph != NULL){
        printf("dt=%f\tNparticles=%d\n", sph->dt, sph->Nparticles);
    }
   
    /* One caveat is that since REBOUNDx does not know the structure definition, custom parameters will not be written or read from REBOUNDx binary files*/

    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    reb_free_simulation(sim);
}
