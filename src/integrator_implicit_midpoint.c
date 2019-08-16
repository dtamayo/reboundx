/**
 * @file    integrator_implicit_midpoint.c
 * @brief   Symplectic numerical integration scheme.
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
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

static void avg_particles(struct reb_particle* const ps_avg, struct reb_particle* const ps1, struct reb_particle* const ps2, int N){
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
}

static int compare(struct reb_particle* ps1, struct reb_particle* ps2, int N){
    double tot2 = 0.;
    double deltatot2 = 0.;
    for(int i=0; i<N; i++){
        const double dvx = ps1[i].vx - ps2[i].vx;
        const double dvy = ps1[i].vy - ps2[i].vy;
        const double dvz = ps1[i].vz - ps2[i].vz;
        deltatot2 += dvx*dvx + dvy*dvy + dvz*dvz;
        tot2 += ps1[i].vx*ps1[i].vx + ps1[i].vy*ps1[i].vy + ps1[i].vz*ps1[i].vz;
    }
    if (deltatot2/tot2 < DBL_EPSILON*DBL_EPSILON){
        return 1;
    }
    else{
        return 0;
    }
}

void rebx_im_free_arrays(struct rebx_extras* rebx, struct rebx_force* force){
    struct reb_particle* const ps_final = rebx_get_param(rebx, force->ap, "ps_final");
    free(ps_final);
    struct reb_particle* const ps_prev = rebx_get_param(rebx, force->ap, "ps_prev");
    free(ps_prev);
    struct reb_particle* const ps_avg = rebx_get_param(rebx, force->ap, "ps_avg");
    free(ps_avg);
}

static struct reb_particle* setup(struct rebx_extras* rebx, struct rebx_force* force, const int N){
    rebx_set_param_pointer(rebx, &force->ap, "free_arrays", rebx_im_free_arrays);
    struct reb_particle* const ps_final = malloc(N*sizeof(*ps_final));
    rebx_set_param_pointer(rebx, &force->ap, "im_ps_final", ps_final);
    struct reb_particle* const ps_prev = malloc(N*sizeof(*ps_prev));
    rebx_set_param_pointer(rebx, &force->ap, "im_ps_prev", ps_prev);
    struct reb_particle* const ps_avg = malloc(N*sizeof(*ps_avg));
    rebx_set_param_pointer(rebx, &force->ap, "im_ps_avg", ps_avg);
    
    return ps_final;
}

void rebx_integrator_implicit_midpoint_integrate(struct reb_simulation* const sim, const double dt, struct rebx_force* const force){
    struct rebx_extras* rebx = sim->extras;
    const int N = sim->N - sim->N_var;
    struct reb_particle* ps_final = rebx_get_param(rebx, force->ap, "im_ps_final");
    if (ps_final == NULL){
        ps_final = setup(rebx, force, N);
    }
    // These should not fail since we check above and setup if not there
    struct reb_particle* const ps_prev = rebx_get_param(rebx, force->ap, "im_ps_prev");
    struct reb_particle* const ps_avg = rebx_get_param(rebx, force->ap, "im_ps_avg");
    struct reb_particle* const ps_orig = sim->particles;
    memcpy(ps_final, sim->particles, N*sizeof(*ps_final));
    memcpy(ps_avg, sim->particles, N*sizeof(*ps_orig));
    int n, converged;
    for(n=0;n<10;n++){
        memcpy(ps_prev, ps_final, N*sizeof(*ps_prev));
        force->update_accelerations(sim, force, ps_avg, N);
        for(int i=0; i<N; i++){
            ps_final[i].vx = ps_orig[i].vx + dt*ps_avg[i].ax;
            ps_final[i].vy = ps_orig[i].vy + dt*ps_avg[i].ay;
            ps_final[i].vz = ps_orig[i].vz + dt*ps_avg[i].az;
        }
        converged = compare(ps_final, ps_prev, N);
        if (converged){
            break;
        }
        avg_particles(ps_avg, ps_orig, ps_final, N);
    }
    const int default_max_iterations = 10;
    if(n==default_max_iterations){
        reb_warning(sim, "REBOUNDx: 10 iterations in integrator_implicit_midpoint.c failed to converge. This is typically because the perturbation is too strong for the current implementation.");
    }
    for(int i=0; i<N; i++){
        sim->particles[i].vx = ps_final[i].vx;
        sim->particles[i].vy = ps_final[i].vy;
        sim->particles[i].vz = ps_final[i].vz;
    }
}
