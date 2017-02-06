/**
 * @file    gr.c
 * @brief   Post-newtonian general relativity corrections arising from a single massive body
 * @author  Pengshuai (Sam) Shi, Dan Tamayo, Hanno Rein <tamayo.daniel@gmail.com>
 *
 * @section     LICENSE
 * Copyright (c) 2015 Pengshuai (Sam) Shi, Dan Tamayo, Hanno Rein
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
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $General Relativity$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 P. Shi, D. Tamayo, H. Rein
 * Implementation Paper    *In progress*
 * Based on                `Anderson et al. 1975 <http://labs.adsabs.harvard.edu/adsabs/abs/1975ApJ...200..221A/>`_.
 * C Example               :ref:`c_example_gr`
 * Python Example          `GeneralRelativity.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GeneralRelativity.ipynb>`_.
 * ======================= ===============================================
 * 
 * This assumes that the masses are dominated by a single central body, and should be good enough for most applications with planets orbiting single stars.
 * It ignores terms that are smaller by of order the mass ratio with the central body.
 * It gets both the mean motion and precession correct, and will be significantly faster than :ref:`gr_full`, particularly with several bodies.
 * Adding this effect to several bodies is NOT equivalent to using gr_full.
 * 
 * **Effect Parameters**
 * 
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * c (double)                   Yes         Speed of light in the units used for the simulation.
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * If no particles have gr_source set, effect will assume the particle at index 0 in the particles array is the source.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * gr_source (int)              No          Flag identifying the particle as the source of perturbations.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_gr(struct reb_simulation* const sim, const double C2, const int no_iteration){
    const int N_real = sim->N - sim->N_var;
    const double G = sim->G;

    struct reb_particle* const ps = malloc(N_real*sizeof(*ps));
    struct reb_particle* const ps_j = malloc(N_real*sizeof(*ps_j));
    memcpy(ps, sim->particles, N_real*sizeof(*ps));
    
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
        double vi2=vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
        const double ri = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        int q = 0;
        double A;
        if (no_iteration){
            A = (0.5*vi2 + 3.*mu/ri)/C2;
        }
        else{
            struct reb_vec3d old_v;
            for(q=0; q<10; q++){
                old_v.x = vi.x;
                old_v.y = vi.y;
                old_v.z = vi.z;
                A = (0.5*vi2 + 3.*mu/ri)/C2;
                vi.x = p.vx/(1.-A);
                vi.y = p.vy/(1.-A);
                vi.z = p.vz/(1.-A);
                vi2 =vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
                const double dvx = vi.x - old_v.x;
                const double dvy = vi.y - old_v.y;
                const double dvz = vi.z - old_v.z;
                if ((dvx*dvx + dvy*dvy + dvz*dvz)/vi2 < DBL_EPSILON*DBL_EPSILON){
                    break;
                }
            }
        }
        //fprintf(stderr, "%d\n", q);
        const double B = (mu/ri - 1.5*vi2)*mu/(ri*ri*ri)/C2;
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
    for (int i=0; i<N_real; i++){
        sim->particles[i].ax += ps[i].ax;
        sim->particles[i].ay += ps[i].ay;
        sim->particles[i].az += ps[i].az;
    }
    free(ps);
    free(ps_j);
    free(eta);
}

void rebx_gr(struct reb_simulation* const sim, struct rebx_effect* const gr){ // First find gr sources
    double* c = rebx_get_param_check(gr, "c", REBX_TYPE_DOUBLE);
    if (c == NULL){
        reb_error(sim, "Need to set speed of light in gr effect.  See examples in documentation.\n");
        return;
    }
    const double C2 = (*c)*(*c);
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    int no_iteration;
    if(rebx_get_param_check(gr, "no_iteration", REBX_TYPE_INT) != NULL){
        no_iteration = 1;
    }
    else{
        no_iteration = 0;
    }
    rebx_calculate_gr(sim, C2, no_iteration);  // gr_source not found, default to index=0
}

static double rebx_calculate_gr_hamiltonian(struct reb_simulation* const sim, const double C2, const int source_index){
    const int N_real = sim->N - sim->N_var;
    const double G = sim->G;

    struct reb_particle* const ps_j = malloc(N_real*sizeof(*ps_j));
    struct reb_particle* const ps = sim->particles; 
    // Calculate Newtonian potentials

    double V_newt = 0.;
    for(int i=0; i<N_real; i++){
        const struct reb_particle pi = ps[i];
        for(int j=i+1; j<N_real; j++){
            const struct reb_particle pj = ps[j];
            const double dx = pi.x - pj.x;
            const double dy = pi.y - pj.y;
            const double dz = pi.z - pj.z;
            const double softening2 = sim->softening*sim->softening;
            const double r2 = dx*dx + dy*dy + dz*dz + softening2;
            V_newt -= G*pi.m*pj.m/sqrt(r2);
        }
    }
   
    // Transform to Jacobi coordinates
    const struct reb_particle source = ps[0];
	const double mu = G*source.m;
    double* const eta = malloc(N_real*sizeof(*eta));
    double* const m_j = malloc(N_real*sizeof(*m_j));
    eta[0] = ps[0].m;
    for (unsigned int i=1;i<N_real;i++){
        eta[i] = eta[i-1] + ps[i].m;
        m_j[i] = ps[i].m*eta[i-1]/eta[i];
    }
    m_j[0] = eta[N_real-1];

    to_jacobi_posvel(ps, ps_j, eta, ps, N_real);

    double T = 0.5*m_j[0]*(ps_j[0].vx*ps_j[0].vx + ps_j[0].vy*ps_j[0].vy + ps_j[0].vz*ps_j[0].vz);
    double V_PN = 0.;
    for (int i=1; i<N_real; i++){
        struct reb_particle p = ps_j[i];
        const double rdoti2 = p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
        double vtildei2 = rdoti2;
        double A, old_vtildei2;
        const double ri = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        const double vscale2 = mu/ri; // characteristic v^2
        for(int q=0; q<10; q++){
            old_vtildei2 = vtildei2;
            A = (0.5*vtildei2 + 3.*vscale2)/C2;
            vtildei2 = rdoti2/((1.-A)*(1.-A));
            if ((vtildei2 - old_vtildei2)/vtildei2 < DBL_EPSILON){
                break;
            }
        }

        V_PN += m_j[i]*(0.5*mu*mu/(ri*ri) - 0.125*vtildei2*vtildei2 - 1.5*mu*vtildei2/ri);
        T += 0.5*m_j[i]*vtildei2;
    }
    V_PN /= C2;
    
    free(ps_j);
    free(eta);
    free(m_j);
    
	return T + V_newt + V_PN;
}

double rebx_gr_hamiltonian(struct reb_simulation* const sim, const struct rebx_effect* const gr){ 
    double* c = rebx_get_param_check(gr, "c", REBX_TYPE_DOUBLE);
    if (c == NULL){
        reb_error(sim, "Need to set speed of light in gr effect.  See examples in documentation.\n");
        return 0;
    }
    const double C2 = (*c)*(*c);
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    for (int i=0; i<N_real; i++){
        if (rebx_get_param_check(&particles[i], "gr_source", REBX_TYPE_INT) != NULL){
            return rebx_calculate_gr_hamiltonian(sim, C2, i);
        }
    }
    return rebx_calculate_gr_hamiltonian(sim, C2, 0);
}

