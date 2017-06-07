/**
 * @file    euler.c
 * @brief   Euler's method
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
#include "rebxtools.h"

void rebx_integrator_euler_integrate(struct reb_simulation* const sim, const double dt, struct rebx_effect* const effect){
    const int N = sim->N - sim->N_var;
    effect->force(sim, effect, sim->particles, N);
    double* Edissipated = rebx_get_param_check(effect, "Edissipated", REBX_TYPE_DOUBLE);
    if(Edissipated != NULL){
        const double Edot = rebx_Edot(sim->particles, N);
        *Edissipated += dt*Edot;
        
    }
    //fprintf(stderr, "%.16e\t%f\n", sim->particles[1].ax, dt);
    //fprintf(stderr, "%.16e\t%f\n", sim->particles[1].ay, dt);
    for(int i=0; i<N; i++){
        sim->particles[i].vx += dt*sim->particles[i].ax;
        sim->particles[i].vy += dt*sim->particles[i].ay;
        sim->particles[i].vz += dt*sim->particles[i].az;
    }
}
