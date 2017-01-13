/**
 * General central force.
 * 
 * This example shows how to add a general central force.
 * If you have GLUT installed for the visualization, press 'w' and/or 'c' for a clearer view of
 * the whole orbit.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    struct rebx_extras* rebx = rebx_init(sim);

    struct reb_particle star = {0};
    star.m     = 1.;   
    reb_add(sim, star);
   
    double m = 1./333000;
    double a = 1.; // put planet close to enhance precession so it's visible in visualization (this would put planet inside the Sun!)
    double e = 0.0167;
    double omega = 0.;
    double f = 0.;
    
    struct reb_particle planet = reb_tools_orbit2d_to_particle(sim->G, star, m, a, e, omega, f);
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    rebx_add(rebx, "moon_quadrupole_laskar");
    /* We first choose a power (must be a double!) for our central force (here F goes as r^-1).
     * We then need to add it to the particle(s) that will act as central sources for this force.*/

    struct reb_particle* ps = sim->particles;

    double* f_mql = rebx_add_param(&ps[1], "f_mql", REBX_TYPE_DOUBLE);
    *f_mql = 0.9473;

/*  old stuff
    double* m_ratio_earthmoon_mql = rebx_add_param(&ps[1], "m_ratio_earthmoon_mql", REBX_TYPE_DOUBLE);
    *m_ratio_earthmoon_mql = 81.3007;

    double* a0_mql = rebx_add_param(&ps[1], "a0_mql", REBX_TYPE_DOUBLE);
    *a0_mql = 0.0025696;

    double* a1_mql = rebx_add_param(&ps[1], "a1_mql", REBX_TYPE_DOUBLE);
    *a1_mql = 0.101773133860118;

    double* a2_mql = rebx_add_param(&ps[1], "a2_mql", REBX_TYPE_DOUBLE);
    *a2_mql = -0.0178457555910623;

    double* alpha_mql = rebx_add_param(&ps[1], "alpha_mql", REBX_TYPE_DOUBLE);
    *alpha_mql = 0.113381622646105; // USE 1/9 AND FIT FOR OTHER PARAMETERS
*/

    /* We can also use the function rebx_central_force_Acentral to calculate the Acentral required
     * for particles[1] (around primary particles[0]) to have a pericenter precession rate of
     * pomegadot, given a gammacentral value: */
    
    double E0 = rebx_moon_quadrupole_laskar_hamiltonian(sim) + reb_tools_energy(sim); // relativistic hamiltonian
    double tmax = 3.e4*2.*M_PI;
    sim->integrator = REB_INTEGRATOR_IAS15;
    reb_integrate(sim, tmax); 
    rebx_free(rebx);    // this explicitly frees all the memory allocated by REBOUNDx 
    double Ef = rebx_moon_quadrupole_laskar_hamiltonian(sim) + reb_tools_energy(sim); // relativistic hamiltonian
    printf("%e\n",(Ef-E0)/E0);
    reb_free_simulation(sim);
}
