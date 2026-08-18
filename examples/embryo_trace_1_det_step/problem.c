#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"


char TITLE[100] = "embryo_trace_";
char TITLE_FAMTREE[100] = "famtree_";
char TITLE_HEARTBEAT[100] = "heartbeat_";
double RHO = 5.05e6; //3 g/cm^3

double get_radii(double m, double rho){
    return pow((3*m)/(4*M_PI*rho),1./3.);
}

void heartbeat(struct reb_simulation* sim){
    if (reb_simulation_output_check(sim, 1)){
        FILE *fp = fopen("heartbeat.txt", "a");  // append mode
        if (fp == NULL){
            perror("Error opening file");
            return;
        }
        struct reb_vec3d angular_momentum = reb_simulation_angular_momentum(sim);
        double total_mass = 0;
        for (int i = 0; i < sim->N; i++){
            total_mass = total_mass + sim->particles[i].m;
        }
        fprintf(fp, "%e,", sim->t);
        fprintf(fp, "%f,", sim->walltime);
        fprintf(fp, "%d,", sim->N);
        fprintf(fp, "%e,", total_mass);
        fprintf(fp, "%e, ", angular_momentum.x);
        fprintf(fp, "%e, ", angular_momentum.y);
        fprintf(fp, "%e, ", angular_momentum.z);
        fprintf(fp, "%e\n", reb_simulation_energy(sim));

        fclose(fp);
    }
    
}

int main(int argc, char* argv[]){
    char* filename = "/home/tajer.1/reboundx/examples/embryo_trace_1_step/embryo_trace_1";
    struct reb_simulationarchive* sa = reb_simulationarchive_create_from_file(filename);
    struct reb_simulation* r = reb_simulation_create_from_simulationarchive(sa,1450);
    reb_simulationarchive_free(sa);
    printf("Found simulationarchive. Loaded snapshot at t=%.16f.\n",r->t);

    //r->integrator = REB_INTEGRATOR_TRACE;
    //r->ri_trace.S_peri = reb_integrator_trace_switch_peri_none;
    //r->ri_trace.r_crit_hill = 3.*1.21;
    //r->collision = REB_COLLISION_DIRECT;
    //r->G = 39.476926421373;
    //r->dt = 6./365.;
    //r->exact_finish_time = 0;
    //r->heartbeat = heartbeat;
    //r->rand_seed = 1;
    //Setting the collision resolve module
    struct rebx_extras* rebx = rebx_attach(r);
    struct rebx_collision_resolve* fragmenting = rebx_load_collision_resolve(rebx, "fragmenting_collisions");
    rebx_add_collision_resolve(rebx, fragmenting);

    //Choose minimum fragment mass
    double min_frag_mass = 1.4e-8;
    rebx_set_param_double(rebx, &fragmenting->ap, "fc_min_frag_mass", min_frag_mass);
    rebx_set_param_pointer(rebx, &fragmenting->ap, "fc_particle_list_file", TITLE_FAMTREE);


    //r->dt = 6./365.;
    //r->heartbeat = heartbeat;

    // Reset auto-save settings inherited from the loaded snapshot
    //r->simulationarchive_auto_interval = 0;
    //r->simulationarchive_auto_walltime = 0;
    //r->simulationarchive_auto_step = 0;
    
    // Now you can safely set the new step-based save
    reb_simulation_save_to_file_step(r, "archive_timestep.bin", 1);
    //reb_simulation_save_to_file_interval(r, "archive_timestep_1.bin", 6./356.);

    reb_simulation_integrate(r, r->t+1000);
    reb_simulation_free(r);

}
