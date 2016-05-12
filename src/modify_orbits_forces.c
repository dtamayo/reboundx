/**
 * @file    modify_orbits_forces.c
 * @brief   Update orbital elements with prescribed timescales using forces.
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "modify_orbits_forces.h"
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

#define TWOPI 6.2831853071795862

struct rebx_params_modify_orbits_forces* rebx_add_modify_orbits_forces(struct rebx_extras* rebx){
	struct rebx_params_modify_orbits_forces* params = malloc(sizeof(*params));
    params->coordinates = JACOBI;
    int force_is_velocity_dependent = 1;
    rebx_add_force(rebx, params, "modify_orbits_forces", rebx_modify_orbits_forces, force_is_velocity_dependent);
    return params;
}

void rebx_modify_orbits_forces(struct reb_simulation* const sim, struct rebx_effect* const effect){
    const struct rebx_params_modify_orbits_forces* const params = effect->paramsPtr;
    const int N_real = sim->N - sim->N_var;

    struct reb_particle com;
    switch(params->coordinates){
    case JACOBI:
        com = reb_get_com(sim);                     // We start with outermost particle, so start with COM and peel off particles
        break;
    case BARYCENTRIC:
        com = reb_get_com(sim);                     // COM of whole system
        break;
    case HELIOCENTRIC:
        com = sim->particles[0];                    // Use the central body as the reference
        break;
    default:
        fprintf(stderr, "coordinates in parameters for modify_orbits_forces are not supported.\n");
        exit(1);
    }

    struct reb_particle* p0 = &sim->particles[0];
    for(int i=N_real-1;i>0;--i){
        struct reb_particle* p = &sim->particles[i];
        if(params->coordinates == JACOBI){
            rebxtools_update_com_without_particle(&com, p);
        }
        double tau_a = rebx_get_param_double(p, "tau_a");
        double tau_e = rebx_get_param_double(p, "tau_e");
        double tau_inc = rebx_get_param_double(p, "tau_inc");

        if(isnan(tau_a)){
            rebx_set_param_double(p, "tau_a", INFINITY);
            tau_a = INFINITY;
        }
        if(isnan(tau_e)){
            rebx_set_param_double(p, "tau_e", INFINITY);
            tau_e = INFINITY;
        }
        if(isnan(tau_inc)){
            rebx_set_param_double(p, "tau_inc", INFINITY);
            tau_inc = INFINITY;
        }

        const double dvx = p->vx - com.vx;
        const double dvy = p->vy - com.vy;
        const double dvz = p->vz - com.vz;

        double ax = 0.;
        double ay = 0.;
        double az = 0.;

        ax =  dvx/(2.*tau_a);
        ay =  dvy/(2.*tau_a);
        az =  dvz/(2.*tau_a);

        if (tau_e < INFINITY || tau_inc < INFINITY){   // need h and e vectors for both types
            const double dx = p->x-com.x;
            const double dy = p->y-com.y;
            const double dz = p->z-com.z;
            const double r = sqrt ( dx*dx + dy*dy + dz*dz );
            const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
            const double prefac = 2*vr/r;
            ax += prefac*dx/tau_e;
            ay += prefac*dy/tau_e;
            az += prefac*dz/tau_e + 2.*dvz/tau_inc;
        }
        p->ax += ax;
        p->ay += ay;
        p->az += az;
    }
    reb_move_to_com(sim);
}

