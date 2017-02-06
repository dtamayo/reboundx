/**
 * @file    core.c
 * @brief   Central internal functions for REBOUNDx (not called by user)
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* Main routines called each timestep. */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "core.h"
#include "rebound.h"
#define STRINGIFY(s) str(s)
#define str(s) #s

const char* rebx_build_str = __DATE__ " " __TIME__; // Date and time build string. 
const char* rebx_version_str = "2.16.0";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* rebx_githash_str = STRINGIFY(GITHASH);             // This line gets updated automatically. Do not edit manually.

/*void rebx_gr_acc(struct rebx_extras* const rebx, double* acc, const double C2){
    struct reb_simulation* const sim = rebx->sim;
    struct reb_particle* const particles = sim->particles;
    const int source_index = 0;
    const int N_real = sim->N - sim->N_var;
    const double G = sim->G;
    const unsigned int _gravity_ignore_10 = sim->gravity_ignore_terms==1;
    
    const double mu = G*particles[source_index].m;
    double aoverm10x, aoverm10y, aoverm10z;
   
    for(int i=0; i<N_real;i++){
        sim->particles[i].ax = 0;
        sim->particles[i].ay = 0;
        sim->particles[i].az = 0;
    }
    if (_gravity_ignore_10){
        const double dx = particles[0].x - particles[1].x;
        const double dy = particles[0].y - particles[1].y;
        const double dz = particles[0].z - particles[1].z;
        const double softening2 = sim->softening*sim->softening;
        const double r2 = dx*dx + dy*dy + dz*dz + softening2;
        const double r = sqrt(r2);
        const double prefac = G/(r2*r);
        
        aoverm10x = prefac*dx;
        aoverm10y = prefac*dy;
        aoverm10z = prefac*dz;
    }
    
    const struct reb_particle source = particles[source_index];
    for (int i=0; i<N_real; i++){
        if(i == source_index){
            continue;
        }
        const struct reb_particle pi = particles[i];
        
        const double dx = pi.x - source.x;
        const double dy = pi.y - source.y;
        const double dz = pi.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
        const double vx = pi.vx;
        const double vy = pi.vy;
        const double vz = pi.vz;
        const double v2 = vx*vx + vy*vy + vz*vz;
        double ax = pi.ax;
        double ay = pi.ay;
        double az = pi.az;
        if(_gravity_ignore_10 && i==1){
            ax += particles[0].m*aoverm10x;
            ay += particles[0].m*aoverm10y;
            az += particles[0].m*aoverm10z;
        }
        if(_gravity_ignore_10 && i==0){
            ax -= particles[1].m*aoverm10x;
            ay -= particles[1].m*aoverm10y;
            az -= particles[1].m*aoverm10z;
        }
        
        const double a1_x = (mu*mu*dx/(r2*r2) - 3.*mu*v2*dx/(2.*r2*r))/C2;
        const double a1_y = (mu*mu*dy/(r2*r2) - 3.*mu*v2*dy/(2.*r2*r))/C2;
        const double a1_z = (mu*mu*dz/(r2*r2) - 3.*mu*v2*dz/(2.*r2*r))/C2;
        
        const double va = vx*ax + vy*ay + vz*az;
        const double rv = dx*vx + dy*vy + dz*vz;
        
        ax = a1_x-(va*vx + v2*ax/2. + 3.*mu*(ax*r-vx*rv/r)/r2)/C2;
        ay = a1_y-(va*vy + v2*ay/2. + 3.*mu*(ay*r-vy*rv/r)/r2)/C2;
        az = a1_z-(va*vz + v2*az/2. + 3.*mu*(az*r-vz*rv/r)/r2)/C2;
        
        particles[i].ax += ax;
        particles[i].ay += ay;
        particles[i].az += az;
        particles[source_index].ax -= pi.m/source.m*ax;
        particles[source_index].ay -= pi.m/source.m*ay; 
        particles[source_index].az -= pi.m/source.m*az; 
    }	

    acc[0] = particles[0].ax;
    acc[1] = particles[0].ay;
    acc[2] = particles[0].az;
    acc[3] = particles[1].ax;
    acc[4] = particles[1].ay;
    acc[5] = particles[1].az;
}
*/

static void gr(struct reb_simulation* const sim, struct reb_particle* ps, const double C2){
    const int N_real = sim->N - sim->N_var;
    const double G = sim->G;
    
    struct reb_particle* const ps_j = malloc(N_real*sizeof(*ps_j));
    
    // Calculate Newtonian accelerations
    for(int i=0; i<N_real; i++){
        ps[i].ax = 0.;
        ps[i].ay = 0.;
        ps[i].az = 0.;
    }
    
    for(int i=0; i<N_real; i++){
        const struct reb_particle pi = ps[i];
        for(int j=i+1; j<N_real; j++){
            const struct reb_particle pj = ps[j];
            const double dx = pi.x - pj.x;
            const double dy = pi.y - pj.y;
            const double dz = pi.z - pj.z;
            const double softening2 = sim->softening*sim->softening;
            const double r2 = dx*dx + dy*dy + dz*dz + softening2;
            const double r = sqrt(r2);
            const double prefac = G/(r2*r);
            ps[i].ax -= prefac*pj.m*dx;
            ps[i].ay -= prefac*pj.m*dy;
            ps[i].az -= prefac*pj.m*dz;
            ps[j].ax += prefac*pi.m*dx;
            ps[j].ay += prefac*pi.m*dy;
            ps[j].az += prefac*pi.m*dz;
        }
    }
    
    // Transform to Jacobi coordinates
    const struct reb_particle source = ps[0];
    const double mu = G*source.m;
    double* const eta = malloc(N_real*sizeof(*eta));
    eta[0] = ps[0].m;
    for (unsigned int i=1;i<N_real;i++){
        eta[i] = eta[i-1] + ps[i].m;
    }
    
    to_jacobi_posvel(ps, ps_j, eta, ps, N_real);
    to_jacobi_acc(ps, ps_j, eta, ps, N_real);
    
    for (int i=1; i<N_real; i++){
        struct reb_particle p = ps_j[i];
        struct reb_vec3d vi;
        vi.x = p.vx;
        vi.y = p.vy;
        vi.z = p.vz;
        double vi2, A;
        const double ri = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        for(int q=0; q<10; q++){
            vi2 = vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
            A = (0.5*vi2 + 3.*mu/ri)/C2;
            vi.x = p.vx/(1.-A);
            vi.y = p.vy/(1.-A);
            vi.z = p.vz/(1.-A);
        }
        
        const double B = (mu/ri - 1.5*vi2)*mu/(ri*ri*ri)/C2;
        const double vdotr = vi.x*p.x + vi.y*p.x + vi.z*p.z;
        const double rdotrdot = p.x*p.vx + p.y*p.vy + p.z*p.vz;
        
        struct reb_vec3d vidot;
        vidot.x = p.ax + B*p.x;
        vidot.y = p.ay + B*p.y;
        vidot.z = p.az + B*p.z;
        
        const double vdota = vi.x*vidot.x + vi.y*vidot.y + vi.z*vidot.z;
        const double D = (vdota - 3.*mu/(ri*ri*ri)*rdotrdot)/C2;
        ps_j[i].ax = B*(1.-A)*p.x - A*p.ax - D*vi.x;
        ps_j[i].ay = B*(1.-A)*p.y - A*p.ay - D*vi.y;
        ps_j[i].az = B*(1.-A)*p.z - A*p.az - D*vi.z;
    }
    ps_j[0].ax = 0.;
    ps_j[0].ay = 0.;
    ps_j[0].az = 0.;
    
    to_inertial_acc(ps, ps_j, eta, ps, N_real);
}

void drag_force(struct reb_simulation* const r, struct reb_particle* const ps){
    double tau = 100.;
    struct reb_particle* source = &ps[0];
    struct reb_particle* p = &ps[1];
    const double dvx = p->vx - source->vx;
    const double dvy = p->vy - source->vy;
    const double dvz = p->vz - source->vz;
    p->ax -= dvx/tau;
    p->ay -= dvy/tau;
    p->az -= dvz/tau;
    source->ax += p->m/source->m*dvx/tau;
    source->ay += p->m/source->m*dvy/tau;
    source->az += p->m/source->m*dvz/tau;
}

void avg_particles(struct reb_particle* const ps_avg, struct reb_particle* const ps1, struct reb_particle* const ps2, int N){
    for(int i=0; i<N; i++){
        ps_avg[i].x = 0.5*(ps1[i].x + ps2[i].x);
        ps_avg[i].y = 0.5*(ps1[i].y + ps2[i].y);
        ps_avg[i].z = 0.5*(ps1[i].z + ps2[i].z); 
        ps_avg[i].vx = 0.5*(ps1[i].vx + ps2[i].vx);
        ps_avg[i].vy = 0.5*(ps1[i].vy + ps2[i].vy);
        ps_avg[i].vz = 0.5*(ps1[i].vz + ps2[i].vz);
        ps_avg[i].ax = 0.;
        ps_avg[i].ay = 0.;
        ps_avg[i].az = 0.;
        ps_avg[i].m = 0.5*(ps1[i].m + ps2[i].m);
    }
    //fprintf(stderr, "ps = %.16f\t ps_old = %.16f\t Avg = %.16f\n", ps1[1].vy, ps2[1].vy, ps_avg[1].vy);
}

int compare(struct reb_particle* ps1, struct reb_particle* ps2, int N){
    //fprintf(stderr, "ps = %.16f\t ps_old = %.16f\t%.4e\n", ps1[1].vy, ps2[1].vy, fabs((ps1[1].vy - ps2[1].vy)/ps1[1].vy));
    for(int i=0; i<N; i++){
        if (ps1[i].vx != ps2[i].vx || ps1[i].vy != ps2[i].vy || ps1[i].vz != ps2[i].vz){
            return 0;
        }
    }
    return 1;
}

void integrate(struct reb_simulation* const r, const double dt){
    const int N = r->N - r->N_var;
    r->ri_whfast.recalculate_jacobi_this_timestep = 1;
    struct reb_particle* const ps = malloc(N*sizeof(*ps));
    memcpy(ps, r->particles, r->N*sizeof(*ps));
    struct reb_particle* ps_orig = malloc(r->N*sizeof(*ps_orig));
    struct reb_particle* ps_old = malloc(r->N*sizeof(*ps_orig));
    struct reb_particle* ps_avg = malloc(r->N*sizeof(*ps_avg));
    memcpy(ps_orig, r->particles, r->N*sizeof(*ps_orig));

    const double C2 = 1.e6;
    int n, converged;
    for(n=1;n<10;n++){
        //memcpy(ps_old, ps, r->N*sizeof(*ps_old));
        avg_particles(ps_avg, ps_orig, ps, N);
        gr(r, ps_avg, C2);//drag_force(r, ps_avg);
        //fprintf(stderr, "ps_orig = %.16f\n", ps_avg[1].ay);
        for(int i=0; i<N; i++){
            ps[i].vx = ps_orig[i].vx + dt*ps_avg[i].ax;
            ps[i].vy = ps_orig[i].vy + dt*ps_avg[i].ay;
            ps[i].vz = ps_orig[i].vz + dt*ps_avg[i].az;
        }
        //converged = compare(ps, ps_old, N);
        //if (converged){
        //    break;
        //}
    }
    //fprintf(stderr, "%d\n", n);
    //fprintf(stderr, "%e\n", fabs((ps[1].vy-ps_orig[1].vy)/ps_orig[1].vy));
    //double v = sqrt(ps[1].vx*ps[1].vx + ps[1].vy*ps[1].vy);
    //double v_orig = sqrt(ps_orig[1].vx*ps_orig[1].vx + ps_orig[1].vy*ps_orig[1].vy);
    //fprintf(stderr, "%e\n", fabs((v-v_orig)/v_orig));
    for(int i=0; i<N; i++){
        r->particles[i].vx = ps[i].vx;
        r->particles[i].vy = ps[i].vy;
        r->particles[i].vz = ps[i].vz;
    }
}

/*****************************
 Initialization routines.
 ****************************/

struct rebx_extras* rebx_init(struct reb_simulation* sim){  // reboundx.h
    struct rebx_extras* rebx = malloc(sizeof(*rebx));
    rebx_initialize(sim, rebx);
    return rebx;
}

void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx){
    sim->extras = rebx;
    sim->ri_whfast.symplectic_operator = integrate;
    rebx->sim = sim;
	rebx->params_to_be_freed = NULL;
	rebx->effects = NULL;
    
    if(sim->additional_forces || sim->post_timestep_modifications){
        reb_warning(sim, "REBOUNDx will overwrite sim->additional_forces and sim->post_timestep_modifications.  If you are using custom functions you wrote and setting them through these function pointers, this won't work.  You should instead add your custom effects through REBOUNDx.  See http://reboundx.readthedocs.org/en/latest/c_examples.html#adding-custom-post-timestep-modifications-and-forces for a tutorial.");
    }
}

/*****************************
 Garbage Collection Routines
 ****************************/

void rebx_remove_from_simulation(struct reb_simulation* sim){
    sim->additional_forces = NULL;
    sim->post_timestep_modifications = NULL;
}

void rebx_free(struct rebx_extras* rebx){                   // reboundx.h
    rebx_free_params(rebx);
    rebx_free_effects(rebx);
    free(rebx);
}

void rebx_free_params(struct rebx_extras* rebx){
    struct rebx_param_to_be_freed* current = rebx->params_to_be_freed;
    struct rebx_param_to_be_freed* temp_next;
    while(current != NULL){
        temp_next = current->next;
        free(current->param->contents);
        free(current->param);
        free(current);
        current = temp_next;
    }
}

void rebx_free_effects(struct rebx_extras* rebx){
    struct rebx_effect* current = rebx->effects;
    struct rebx_effect* temp_next;

    while(current != NULL){
        temp_next = current->next;
        free(current);
        current = temp_next;
    }
}

/**********************************************
 Functions executing forces & ptm each timestep
 *********************************************/

void rebx_forces(struct reb_simulation* sim){
	struct rebx_extras* rebx = sim->extras;
	struct rebx_effect* current = rebx->effects;
	while(current != NULL){
        if(current->force != NULL){
		    current->force(sim, current);
        }
		current = current->next;
	}
}

void rebx_post_timestep_modifications(struct reb_simulation* sim){
	struct rebx_extras* rebx = sim->extras;
	struct rebx_effect* current = rebx->effects;
	while(current != NULL){
        if(current->ptm != NULL){
		    current->ptm(sim, current);
        }
		current = current->next;
	}
}

/**********************************************
 Adders for linked lists in extras structure
 *********************************************/

static struct rebx_effect* rebx_add_effect(struct rebx_extras* rebx, const char* name){
    struct rebx_effect* effect = malloc(sizeof(*effect));
    
    effect->hash = reb_hash(name);
    effect->ap = NULL;
    effect->force = NULL;
    effect->ptm = NULL;
    effect->rebx = rebx;
    
    effect->name = malloc(strlen(name) + 1); // +1 for \0 at end
    if (effect->name != NULL){
        strcpy(effect->name, name);
    }
    
    effect->next = rebx->effects;
    rebx->effects = effect;
   
    struct reb_particle* p = (struct reb_particle*)effect;
    p->ap = (void*)REBX_OBJECT_TYPE_EFFECT;
    return effect;
}

struct rebx_effect* rebx_add(struct rebx_extras* rebx, const char* name){
    struct rebx_effect* effect = rebx_add_effect(rebx, name);
    struct reb_simulation* sim = rebx->sim;

    if (effect->hash == reb_hash("modify_orbits_direct")){
        effect->ptm = rebx_modify_orbits_direct;
    }
    else if (effect->hash == reb_hash("modify_orbits_forces")){
        sim->force_is_velocity_dependent = 1;
        effect->force = rebx_modify_orbits_forces;
    }
    else if(effect->hash == reb_hash("gr")){
        sim->force_is_velocity_dependent = 1;
        effect->force = rebx_gr;
    }
    else if (effect->hash == reb_hash("gr_full")){
        sim->force_is_velocity_dependent = 1;
        effect->force = rebx_gr_full;
    }
    else if (effect->hash == reb_hash("gr_potential")){
        effect->force = rebx_gr_potential;
    }
    else if (effect->hash == reb_hash("modify_mass")){
        effect->ptm = rebx_modify_mass;
    }
    else if (effect->hash == reb_hash("radiation_forces")){
        sim->force_is_velocity_dependent = 1;
        effect->force = rebx_radiation_forces;
    }
    else if (effect->hash == reb_hash("tides_precession")){
        effect->force = rebx_tides_precession;
    }
    else if (effect->hash == reb_hash("central_force")){
        effect->force = rebx_central_force;
    }
    else if (effect->hash == reb_hash("track_min_distance")){
        effect->ptm = rebx_track_min_distance;
    }
    else{
        char str[100]; 
        sprintf(str, "Effect '%s' passed to rebx_add not found.\n", name);
        reb_error(sim, str);
    }
    
    // Make sure function pointers in simulation are set, but only set them if there's at least one REBOUNDx effect using them.
    if(effect->force){
        sim->additional_forces = rebx_forces;
    }
    if(effect->ptm){
        sim->post_timestep_modifications = rebx_post_timestep_modifications;
    }
    
    return effect;
}

struct rebx_effect* rebx_add_custom_force(struct rebx_extras* rebx, const char* name, void (*custom_force)(struct reb_simulation* const sim, struct rebx_effect* const effect), const int force_is_velocity_dependent){
    struct rebx_effect* effect = rebx_add_effect(rebx, name);
    effect->force = custom_force;
    struct reb_simulation* sim = rebx->sim;
    if(force_is_velocity_dependent){
        sim->force_is_velocity_dependent = 1;
    }
    rebx->sim->additional_forces = rebx_forces;
    return effect;
}

struct rebx_effect* rebx_add_custom_post_timestep_modification(struct rebx_extras* rebx, const char* name, void (*custom_ptm)(struct reb_simulation* const sim, struct rebx_effect* const effect)){
    struct rebx_effect* effect = rebx_add_effect(rebx, name);
    effect->ptm = custom_ptm;
    rebx->sim->post_timestep_modifications = rebx_post_timestep_modifications;
    return effect;
}
    
void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param){
    struct rebx_param_to_be_freed* newparam = malloc(sizeof(*newparam));
    newparam->param = param;

    newparam->next = rebx->params_to_be_freed;
    rebx->params_to_be_freed = newparam;
}

/*********************************************************************************
 Internal functions for dealing with parameters
 ********************************************************************************/

static enum rebx_object_type rebx_get_object_type(const void* const object){
    struct reb_particle* p = (struct reb_particle*)object;
    if (p->ap == (void*)(intptr_t)REBX_OBJECT_TYPE_EFFECT){    // In add_effect we cast effect to a particle and set p->ap to REBX_OBJECT_TYPE_EFFECT
        return REBX_OBJECT_TYPE_EFFECT;
    }
    else{
        return REBX_OBJECT_TYPE_PARTICLE;
    }
}


// get simulation pointer from an object
static struct reb_simulation* rebx_get_sim(const void* const object){
    switch(rebx_get_object_type(object)){
        case REBX_OBJECT_TYPE_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            return effect->rebx->sim;
        }
        case REBX_OBJECT_TYPE_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            return p->sim;
        }
    }
    return NULL;
} 

struct rebx_param* rebx_create_param(){
    struct rebx_param* newparam = malloc(sizeof(*newparam));
    int collision=0;
    for(int j=REBX_OBJECT_TYPE_EFFECT; j<=REBX_OBJECT_TYPE_PARTICLE; j++){
        /* need this cast to avoid warnings (we are converting 32 bit int enum to 64 bit pointer). I think behavior
         * is implemetation defined, but this is OK because we are always checking the ap pointer in this way, so 
         * just need that whatever implementation defined result we get from this comparison is false for the address
         * of newparam, which will be the head node in the linked list at object->ap*/
        if(newparam == (void*)(intptr_t)j){ 
            collision=1;
        }
    }
    
    if (collision){
        free(newparam);
        struct rebx_param* address[5] = {NULL};
        int i;
        for(i=0; i<=5; i++){
            address[i] = malloc(sizeof(*newparam));
            collision = 0;
            for(int j=REBX_OBJECT_TYPE_EFFECT; j<=REBX_OBJECT_TYPE_PARTICLE; j++){
                if(address[i] == (void*)(intptr_t)j){
                    collision=1;
                }
            }
            if (collision == 0){
                newparam = address[i];
                break;
            }
        }

        if (i==5){
            fprintf(stderr, "REBOUNDx Error:  Can't allocate valid memory for parameter.\n");
            return NULL;
        }
        for(int j=0;j<i;j++){
            free(address[j]);
        }
    }
    return newparam;
}

/*********************************************************************************
 User interface for parameters
 ********************************************************************************/
int rebx_remove_param(const void* const object, const char* const param_name){
    // TODO free memory for deleted node
    uint32_t hash = reb_hash(param_name);
    struct rebx_param* current;
    switch(rebx_get_object_type(object)){
        case REBX_OBJECT_TYPE_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            current = effect->ap;
            if(current->hash == hash){
                effect->ap = current->next;
                return 1;
            }
            break;
        }
        case REBX_OBJECT_TYPE_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            current = p->ap;
            if(current->hash == hash){
                p->ap = current->next;
                return 1;
            }
            break;
        }
    }    
  
    while(current->next != NULL){
        if(current->next->hash == hash){
            current->next = current->next->next;
            return 1;
        }
        current = current->next;
    }
    return 0;
}

size_t rebx_sizeof(enum rebx_param_type param_type){
    switch(param_type){
        case REBX_TYPE_DOUBLE:
        {
            return sizeof(double);
        }
        case REBX_TYPE_INT:
        {
            return sizeof(int);
        }
        case REBX_TYPE_UINT32:
        {
            return sizeof(uint32_t);
        }
        case REBX_TYPE_ORBIT:
        {
            return sizeof(struct reb_orbit);
        }
        case REBX_TYPE_LONGLONG:
        {
            return sizeof(long long);
        }
    }
    return 0; // type not found
}

struct rebx_param* rebx_attach_param_node(void* const object, struct rebx_param* param){
    void* ptr = rebx_get_param(object, param->name);
    if (ptr != NULL){
        char str[300];
        sprintf(str, "REBOUNDx Error: Parameter '%s' passed to rebx_add_param already exists.\n", param->name);
        reb_error(rebx_get_sim(object), str);
        return NULL;
    }
    
    switch(rebx_get_object_type(object)){
        case REBX_OBJECT_TYPE_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            param->next = effect->ap;
            effect->ap = param;
            break;
        }
        case REBX_OBJECT_TYPE_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            param->next = p->ap;
            p->ap = param;
            break;
        }
    }
    
    return param;
}

struct rebx_param* rebx_add_param_node(void* const object, const char* const param_name, enum rebx_param_type param_type, const int ndim, const int* const shape){
    struct rebx_param* newparam = rebx_create_param();
    newparam->name = malloc(strlen(param_name) + 1); // +1 for \0 at end
    if (newparam->name != NULL){
        strcpy(newparam->name, param_name);
    }
    newparam->hash = reb_hash(param_name);
    newparam->param_type = param_type;
    newparam->python_type = -1; // not used by C
    newparam->ndim = ndim;
    newparam->shape = NULL;
    newparam->size = 1;
    if (ndim > 0){
	    size_t shapesize = sizeof(int)*ndim;
	    newparam->shape = malloc(shapesize);
        newparam->strides = malloc(shapesize);
	    memcpy(newparam->shape, shape, shapesize);
	    for(int i=ndim-1;i>=0;i--){ // going backward allows us to calculate strides at the same time
            newparam->strides[i] = newparam->size;  // stride[i] is equal to the product of the shapes for all indices > i
		    newparam->size *= shape[i];
	    }
    }
    size_t element_size = rebx_sizeof(param_type);
    if (element_size){
        newparam->contents = malloc(element_size*newparam->size); // newparam->size = number of elements in array (1 if scalar)
    }

    else{
        char str[300];
        sprintf(str, "REBOUNDx Error: Parameter type '%d' passed to rebx_add_param_node not supported.\n", param_type);
        reb_error(rebx_get_sim(object), str);
    }

    newparam = rebx_attach_param_node(object, newparam);
    return newparam;
}

void* rebx_add_param_array(void* const object, const char* const param_name, enum rebx_param_type param_type, const int ndim, const int* const shape){
    return rebx_add_param_node(object, param_name, param_type, ndim, shape)->contents;
}

void* rebx_add_param(void* const object, const char* const param_name, enum rebx_param_type param_type){
	int ndim=0;
	return rebx_add_param_array(object, param_name, param_type, ndim, NULL);
}

void* rebx_get_param(const void* const object, const char* const param_name){
    struct rebx_param* node = rebx_get_param_node(object, param_name);
    if (node == NULL){
        return NULL;
    }
    return node->contents;
}


struct rebx_param* rebx_get_param_node(const void* const object, const char* const param_name){
    struct rebx_param* current;
    switch(rebx_get_object_type(object)){
        case REBX_OBJECT_TYPE_EFFECT:
        {
            struct rebx_effect* effect = (struct rebx_effect*)object;
            current = effect->ap;
            break;
        }
        case REBX_OBJECT_TYPE_PARTICLE:
        {
            struct reb_particle* p = (struct reb_particle*)object;
            current = p->ap;
            break;
        }
    }
    
    uint32_t hash = reb_hash(param_name);
    while(current != NULL){
        if(current->hash == hash){
            return current; 
        }
        current = current->next;
    }
   
    if (current == NULL){   // param_name not found.  Return immediately.
        return NULL;
    }

    return current;
}

void* rebx_get_param_check(const void* const object, const char* const param_name, enum rebx_param_type param_type){
    struct rebx_param* node = rebx_get_param_node(object, param_name);
    if (node == NULL){
        return NULL;
    }
    
    if (node->param_type != param_type){
        char str[300];
        sprintf(str, "REBOUNDx Error: Parameter '%s' passed to rebx_get_param_check was found but was of wrong type.  See documentation for your particular effect.  In python, you might need to add a dot at the end of the number when assigning a parameter that REBOUNDx expects as a float.\n", param_name);
        reb_error(rebx_get_sim(object), str);
        return NULL;
    }

    return node->contents;
}

struct rebx_effect* rebx_get_effect(struct rebx_extras* const rebx, const char* const effect_name){
    struct rebx_effect* current = rebx->effects;
    uint32_t hash = reb_hash(effect_name);
    while(current != NULL){
        if(current->hash == hash){
            return current;
        }
        current = current->next;
    }
    
    if (current == NULL){   // effect_name not found.  Return immediately.
        return NULL;
    }
    
    return current;
}

/***********************************************************************************
 * Miscellaneous Functions
***********************************************************************************/
double install_test(void){
    struct reb_simulation* sim = reb_create_simulation();
    struct reb_particle p = {0};
    p.m = 1.; 
    reb_add(sim, p); 
    struct reb_particle p1 = reb_tools_orbit2d_to_particle(sim->G, p, 0., 1., 0.2, 0., 0.);
    reb_add(sim, p1);
    reb_integrate(sim, 1.);
    return sim->particles[1].x;
}
