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
#include "core.h"

void rebx_rk4_free_arrays(struct rebx_extras* rebx, struct rebx_force* force){
    struct reb_particle* const k2 = rebx_get_param(rebx, force->ap, "rk4_k2");
    free(k2);
    struct reb_particle* const k3 = rebx_get_param(rebx, force->ap, "rk4_k3");
    free(k3);
}

void rebx_integrator_rk4_integrate(struct reb_simulation* const sim, const double dt, struct rebx_force* const force){
    struct rebx_extras* rebx = sim->extras;
    const int N = sim->N - sim->N_var;
    rebx_reset_accelerations(sim->particles, N);
    struct reb_particle* k2 = rebx_get_param(rebx, force->ap, "rk4_k2");
    if (k2 == NULL){
        k2 = malloc(N*sizeof(*k2));
        struct reb_particle* k3 = malloc(N*sizeof(*k3));
        rebx_set_param_pointer(rebx, &force->ap, "rk4_k2", k2);
        rebx_set_param_pointer(rebx, &force->ap, "rk4_k3", k3);
        rebx_set_param_pointer(rebx, &force->ap, "free_arrays", rebx_rk4_free_arrays);
        
    }
    struct reb_particle* k3 = rebx_get_param(rebx, force->ap, "rk4_k3");
    memcpy(k2, sim->particles, N*sizeof(*k2));
    memcpy(k3, sim->particles, N*sizeof(*k3));
    
    const double dt2 = dt/2.;
    force->update_accelerations(sim, force, sim->particles, N);  // k1 = sim.particles.a
    
    for(int i=0; i<N; i++){
        k2[i].vx = sim->particles[i].vx + dt2*sim->particles[i].ax;
        k2[i].vy = sim->particles[i].vy + dt2*sim->particles[i].ay;
        k2[i].vz = sim->particles[i].vz + dt2*sim->particles[i].az;
    }
    force->update_accelerations(sim, force, k2, N);
    
    for(int i=0; i<N; i++){
        k3[i].vx = sim->particles[i].vx + dt2*k2[i].ax;
        k3[i].vy = sim->particles[i].vy + dt2*k2[i].ay;
        k3[i].vz = sim->particles[i].vz + dt2*k2[i].az;
    }
    force->update_accelerations(sim, force, k3, N);
    
    for(int i=0; i<N; i++){     // store k2+k3 in k3 and reuse k2 for k4 to avoid a memcpy
        k2[i].vx = sim->particles[i].vx + dt*k3[i].ax;
        k2[i].vy = sim->particles[i].vy + dt*k3[i].ay;
        k2[i].vz = sim->particles[i].vz + dt*k3[i].az;
        k3[i].ax += k2[i].ax;
        k3[i].ay += k2[i].ay;
        k3[i].az += k2[i].az;
    }
    rebx_reset_accelerations(k2, N);
    force->update_accelerations(sim, force, k2, N);
    
    const double dt6 = dt/6.;
    for(int i=0; i<N; i++){
        sim->particles[i].vx += dt6*(sim->particles[i].ax + k2[i].ax + 2.*k3[i].ax);
        sim->particles[i].vy += dt6*(sim->particles[i].ay + k2[i].ay + 2.*k3[i].ay);
        sim->particles[i].vz += dt6*(sim->particles[i].az + k2[i].az + 2.*k3[i].az);
    }
}
