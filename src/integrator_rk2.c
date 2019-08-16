/**
 * @file    rk2.c
 * @brief   2nd order Runge Kutta method (Ralston's method)
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

void rebx_rk2_free_arrays(struct rebx_extras* rebx, struct rebx_force* force){
    struct reb_particle* const k2 = rebx_get_param(rebx, force->ap, "rk2_k2");
    free(k2);
}

void rebx_integrator_rk2_integrate(struct reb_simulation* const sim, const double dt, struct rebx_force* const force){
    struct rebx_extras* rebx = sim->extras;
    const int N = sim->N - sim->N_var;
    struct reb_particle* k2 = rebx_get_param(rebx, force->ap, "rk2_k2");
    if (k2 == NULL){
        k2 = malloc(N*sizeof(*k2));
        rebx_set_param_pointer(rebx, &force->ap, "rk2_k2", k2);
        rebx_set_param_pointer(rebx, &force->ap, "free_arrays", rebx_rk2_free_arrays);
        
    }
    memcpy(k2, sim->particles, N*sizeof(*k2));

    force->update_accelerations(sim, force, sim->particles, N);
    const double a21 = 2.*dt/3.;
    for(int i=0; i<N; i++){
        k2[i].vx = sim->particles[i].vx + a21*sim->particles[i].ax;
        k2[i].vy = sim->particles[i].vy + a21*sim->particles[i].ay;
        k2[i].vz = sim->particles[i].vz + a21*sim->particles[i].az;
    }

    force->update_accelerations(sim, force, k2, N);

    const double b1 = dt/4.;
    const double b2 = 3.*dt/4.;
    for(int i=0; i<N; i++){
        sim->particles[i].vx += b1*sim->particles[i].ax + b2*k2[i].ax;
        sim->particles[i].vy += b1*sim->particles[i].ay + b2*k2[i].ay;
        sim->particles[i].vz += b1*sim->particles[i].az + b2*k2[i].az;
    }
}
