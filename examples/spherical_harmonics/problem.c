/**
 * Adding a General Spherical Harmonics Model (C22, S22) to a central body.
 *
 * This example demonstrates how to set up the general spherical harmonics
 * model in REBOUNDx to simulate the longitude libration of a GEO satellite 
 * due to Earth's tesseral harmonics.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "spherical_harmonics.h"


int main(int argc, char* argv[]) {
    // 1. Create the REBOUND simulation
    struct reb_simulation* sim = reb_simulation_create();
    
    // Set units implicitly by providing G and values in SI
    sim->G = 6.67430e-11;       // m^3 / (kg s^2)
    sim->integrator = REB_INTEGRATOR_IAS15;

    // 2. Constants for Earth and GEO
    double m_earth = 5.9722e24; // kg
    double r_earth = 6371e3;    // m
    double R_eq    = 6378136.3; // m (Reference equatorial radius)
    double omega_earth = 7.292115e-5; // rad/s
    //double theta0  = 0.0;       // initial angle
    
    // 3. Add Earth
    struct reb_particle earth = {0};
    earth.m = m_earth;
    earth.r = r_earth;
    earth.hash = reb_hash("Earth");
    reb_simulation_add(sim, earth);

    // 4. Calculate GEO radius and velocity
    double r_geo = pow(sim->G * m_earth / pow(omega_earth, 2), 1.0 / 3.0);
    double lon0 = 0.0; // Radians
    double v_circ = sqrt(sim->G * m_earth / r_geo);
    
    // Add GEO satellite
    struct reb_particle sat = {0};
    sat.m = 0.0;
    sat.x = r_geo * cos(lon0);
    sat.y = r_geo * sin(lon0);
    sat.z = 0.0;
    sat.vx = -v_circ * sin(lon0);
    sat.vy =  v_circ * cos(lon0);
    sat.vz = 0.0;
    sat.hash = reb_hash("Sat");
    reb_simulation_add(sim, sat);
    
    reb_simulation_move_to_com(sim);

    // 5. Attach REBOUNDx and get the force
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* gh = rebx_load_force(rebx, "gravitational_harmonics");
    
    // Set Earth's rotation vector (required for tesseral harmonics to rotate)
    struct reb_vec3d Omega = {
    .x = 0.0,
    .y = 0.0,
    .z = omega_earth
    };
    
    
    rebx_set_param_vec3d(
        rebx,
        (struct rebx_node**)&sim->particles[0].ap,
        "Omega",
        Omega
    );
    // 6. Setup the General Spherical Harmonics Model
    int N = 2; // Maximum degree
    int size = (N + 1) * (N + 2) / 2;
    
    // Allocate memory for C and S coefficients using calloc (initializes to zero)
    double* C = (double*)calloc(size, sizeof(double));
    double* S = (double*)calloc(size, sizeof(double));
    
    // Insert the unnormalized C22 and S22 coefficients 
    // The idx(n,m) function converts the two inputs into the correct index for the dense triangular array.
    C[idx(2, 2)] = 1.574615e-6;
    S[idx(2, 2)] = -0.903872e-6;
    
    // Create the model structure. 
    // Notice that in C, we pass pointers to the C and S arrays.
    // It's necessary to pass the lenght of C ans S: N
    rebx_spherical_harmonics_model* model = rebx_spherical_harmonics_create(N, C, S, sim->G * m_earth, R_eq, 0); // 0 indicates that the coefficients are not normalized
    
    // Assign the model to the central particle (Earth)
    rebx_set_param_pointer(rebx, (struct rebx_node**)&sim->particles[0].ap, "spherical_harmonics_model", model);

    // Add the force to the simulation
    rebx_add_force(rebx, gh);

    // 7. Integrate
    double sidereal_day = 2.0 * M_PI / omega_earth;
    double tmax = sidereal_day * 100.0; // Integrate for 100 days as a test
    
    printf("Integrating up to %.2f sidereal days...\\n", tmax / sidereal_day);
    reb_simulation_integrate(sim, tmax); 
    printf("Integration complete.\\n");

    // 8. Clean up memory
    free(C);
    free(S);
    rebx_spherical_harmonics_free(model); // Crucial step in C to prevent memory leaks!
    rebx_free(rebx);
    reb_simulation_free(sim);
    
    return 0;
}