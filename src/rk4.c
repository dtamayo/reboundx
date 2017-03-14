/**
 * @file    rk4.c
 * @brief   4th order Runge Kutta method
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>, Hanno Rein
 *
 * @section LICENSE
 * Copyright (c) 2017 Dan Tamayo, Hanno Rein
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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "rebound.h"
#include "reboundx.h"
#include "rk4.h"

static void rebx_reset_acc(struct reb_particle* const particles, const int N){
    for(int i=0;i<N;i++){
        particles[i].ax = 0.;
        particles[i].ay = 0.;
        particles[i].az = 0.;
    }
}

static void rebx_rk_k(struct reb_particle* const ps_orig, struct reb_particle* const k){

}

void rebx_integrator_rk4_integrate(struct reb_simulation* const sim, const double dt, struct rebx_effect* const effect){
    
    const int N = sim->N - sim->N_var;
    rebx_reset_acc(sim->particles, N);
    struct reb_particle* const k2 = malloc(N*sizeof(*k2));
    struct reb_particle* const k3 = malloc(N*sizeof(*k3));
    memcpy(k2, sim->particles, N*sizeof(*k2));
    memcpy(k3, sim->particles, N*sizeof(*k3));
    /*double E0 = 0.;
    double* Edissipated = rebx_get_param_check(effect, "Edissipated", REBX_TYPE_DOUBLE);
    if(Edissipated != NULL){
        E0 = reb_tools_energy(sim);
    }
    const double dt2 = dt/2.;
    effect->force(sim, effect, sim->particles, N);  // k1 = sim.particles.a
    
    for(int i=0; i<N; i++){
        k2[i].vx = sim->particles[i].vx + dt2*sim->particles[i].ax;
        k2[i].vy = sim->particles[i].vy + dt2*sim->particles[i].ay;
        k2[i].vz = sim->particles[i].vz + dt2*sim->particles[i].az;
    }
    effect->force(sim, effect, k2, N);
    for(int i=0; i<N; i++){
        k3[i].vx = sim->particles[i].vx + dt2*k2[i].ax;
        k3[i].vy = sim->particles[i].vy + dt2*k2[i].ay;
        k3[i].vz = sim->particles[i].vz + dt2*k2[i].az;
    }
    effect->force(sim, effect, k3, N);
    for(int i=0; i<N; i++){     // store k2+k3 in k3 and reuse k2 for k4 to avoid a memcpy
        k2[i].vx = sim->particles[i].vx + dt*k3[i].ax;
        k2[i].vy = sim->particles[i].vy + dt*k3[i].ay;
        k2[i].vz = sim->particles[i].vz + dt*k3[i].az;
        k3[i].ax += k2[i].ax;
        k3[i].ay += k2[i].ay;
        k3[i].az += k2[i].az;
    }
    rebx_reset_acc(k2, N);
    effect->force(sim, effect, k2, N);
    const double dt6 = dt/6.;
    for(int i=0; i<N; i++){
        sim->particles[i].vx += dt6*(sim->particles[i].ax + k2[i].ax + 2.*k3[i].ax);
        sim->particles[i].vy += dt6*(sim->particles[i].ay + k2[i].ay + 2.*k3[i].ay);
        sim->particles[i].vz += dt6*(sim->particles[i].az + k2[i].az + 2.*k3[i].az);
    }
    if(Edissipated != NULL){
        const double Ef = reb_tools_energy(sim);
        *Edissipated += Ef-E0;
    }*/
    double* const Edissipated = rebx_get_param_check(effect, "Edissipated",REBX_TYPE_DOUBLE);
    
    const double dt2 = dt/2.;
    effect->force(sim, effect, sim->particles, N);  // k1 = sim.particles.a
    
    double Edot1, Edot2, Edot3, Edot4;
    if (Edissipated != NULL){
        Edot1 = rebx_Edot(sim->particles, N);
    }
    for(int i=0; i<N; i++){
        k2[i].vx = sim->particles[i].vx + dt2*sim->particles[i].ax;
        k2[i].vy = sim->particles[i].vy + dt2*sim->particles[i].ay;
        k2[i].vz = sim->particles[i].vz + dt2*sim->particles[i].az;
    }
    effect->force(sim, effect, k2, N);
    if (Edissipated != NULL){
        Edot2 = rebx_Edot(k2, N);
    }
    for(int i=0; i<N; i++){
        k3[i].vx = sim->particles[i].vx + dt2*k2[i].ax;
        k3[i].vy = sim->particles[i].vy + dt2*k2[i].ay;
        k3[i].vz = sim->particles[i].vz + dt2*k2[i].az;
    }
    effect->force(sim, effect, k3, N);
    if (Edissipated != NULL){
        Edot3 = rebx_Edot(k3, N);
    }
    for(int i=0; i<N; i++){     // store k2+k3 in k3 and reuse k2 for k4 to avoid a memcpy
        k2[i].vx = sim->particles[i].vx + dt*k3[i].ax;
        k2[i].vy = sim->particles[i].vy + dt*k3[i].ay;
        k2[i].vz = sim->particles[i].vz + dt*k3[i].az;
        k3[i].ax += k2[i].ax;
        k3[i].ay += k2[i].ay;
        k3[i].az += k2[i].az;
    }
    rebx_reset_acc(k2, N);
    effect->force(sim, effect, k2, N);
    if (Edissipated != NULL){
        Edot4 = rebx_Edot(k2, N);
    }
    const double dt6 = dt/6.;
    for(int i=0; i<N; i++){
        sim->particles[i].vx += dt6*(sim->particles[i].ax + k2[i].ax + 2.*k3[i].ax);
        sim->particles[i].vy += dt6*(sim->particles[i].ay + k2[i].ay + 2.*k3[i].ay);
        sim->particles[i].vz += dt6*(sim->particles[i].az + k2[i].az + 2.*k3[i].az);
    }
    if(Edissipated != NULL){
        *Edissipated += dt6*(Edot1 + 2.*(Edot2 + Edot3) + Edot4);
    }
    free(k2);
    free(k3);
}