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
#include "rebound.h"
#include "reboundx.h"
#include "integrator_implicit_midpoint.h"

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
    //fprintf(stderr, "ps = %.16f\t ps_old = %.16f\t Avg = %.16f\n", ps1[1].vy, ps2[1].vy, ps_avg[1].vy);
}

static int compare(struct reb_particle* ps1, struct reb_particle* ps2, int N){
    //fprintf(stderr, "ps = %.16f\t ps_old = %.16f\t%.4e\n", ps1[1].vy, ps2[1].vy, fabs((ps1[1].vy - ps2[1].vy)/ps1[1].vy));
    for(int i=0; i<N; i++){
        if (ps1[i].vx != ps2[i].vx || ps1[i].vy != ps2[i].vy || ps1[i].vz != ps2[i].vz){
            return 0;
        }
    }
    return 1;
}

void rebx_integrator_implicit_midpoint_integrate(struct reb_simulation* const sim, const double dt, struct rebx_effect* const effect){
    const int N = sim->N - sim->N_var;
    struct reb_particle* const ps = malloc(N*sizeof(*ps));
    memcpy(ps, sim->particles, N*sizeof(*ps));
    struct reb_particle* ps_orig = malloc(N*sizeof(*ps_orig));
    struct reb_particle* ps_old = malloc(N*sizeof(*ps_orig));
    struct reb_particle* ps_avg = malloc(N*sizeof(*ps_avg));
    memcpy(ps_orig, sim->particles, N*sizeof(*ps_orig));
    
    int n, converged;
    for(n=1;n<10;n++){
        //memcpy(ps_old, ps, r->N*sizeof(*ps_old));
        avg_particles(ps_avg, ps_orig, ps, N);
        effect->force(sim, effect, ps_avg, N);
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
        sim->particles[i].vx = ps[i].vx;
        sim->particles[i].vy = ps[i].vy;
        sim->particles[i].vz = ps[i].vz;
    }
    free(ps);
    free(ps_orig);
    free(ps_old);
    free(ps_avg);
}
