#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

double gr_hamiltonian(struct reb_simulation* const sim){
    const double C2 = 1.e6;//10064.915*10064.915;
    const source_index = 0;
    const int N_real = sim->N - sim->N_var;
    const double G = sim->G;
    struct reb_particle* const particles = sim->particles;
    const double mu = G*particles[source_index].m;
    
    double e_kin = 0.;
    double e_pot = 0.;
    double e_pn  = 0.;
    const struct reb_particle source = particles[source_index];
    for (int i=0;i<N_real;i++){
        struct reb_particle pi = particles[i];
        if (i != source_index){
            double dx = pi.x - source.x;
            double dy = pi.y - source.y;
            double dz = pi.z - source.z;
            double r2 = dx*dx + dy*dy + dz*dz;
            double r = sqrt(r2);
            
            double vx = pi.vx;
            double vy = pi.vy;
            double vz = pi.vz;
            double v2 = vx*vx + vy*vy + vz*vz;
            
            double A = 1. - (v2/2. + 3.*mu/r)/C2;
            double v_tilde2 = v2/(A*A);
            
            e_kin += 0.5*pi.m*v_tilde2;
            e_pn += (mu*mu*pi.m/(2.*r2) - v_tilde2*v_tilde2*pi.m/8. - 3.*mu*v_tilde2*pi.m/(2.*r))/C2;
        }
        else{
            double source_v2 = source.vx*source.vx + source.vy*source.vy + source.vz*source.vz;
            e_kin += 0.5 * source.m * source_v2;
        }
        for (int j=i+1; j<N_real; j++){
            struct reb_particle pj = particles[j];
            double dx = pi.x - pj.x;
            double dy = pi.y - pj.y;
            double dz = pi.z - pj.z;	
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            
            e_pot -= G*pi.m*pj.m/r;
        }
    }
    return e_kin + e_pot + e_pn;
}

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    struct rebx_extras* rebx = rebx_init(sim);

    struct reb_particle star = {0};
    star.m     = 1.;   
    star.hash  = reb_hash("star");
    reb_add(sim, star);

    double m = 1.e-3;
    double a = 1.;
    double e = 0.1;
    double omega = 0.;
    double f = 0.5;
    
    struct reb_particle planet = reb_tools_orbit2d_to_particle(sim->G, star, m, a, e, omega, f);
    planet.hash = reb_hash("planet");
    reb_add(sim, planet);
    reb_move_to_com(sim);
    
    sim->integrator     = REB_INTEGRATOR_WHFAST;
    sim->dt = 1.e-3;
    double E0 = gr_hamiltonian(sim);
    reb_step(sim);
    //reb_integrate(sim, 10.);
    double Ef = gr_hamiltonian(sim);
    printf("Error: %e\n", fabs((Ef-E0)/E0));
}
