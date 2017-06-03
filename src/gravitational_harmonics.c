/** * @file gravitational_harmonics.c
 * @brief   Add J2n gravitational harmonics to particles
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
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Gravity Fields$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo
 * Implementation Paper    *In progress*
 * Based on                None
 * C Example               :ref:`c_example_J2`
 * Python Example          `J2.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/J2.ipynb>`_.
 * ======================= ===============================================
 * 
 *
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * J2 (double)                  No          J2 coefficient
 * J4 (double)                  No          J4 coefficient
 * R_eq (double)                Yes         Equatorial radius of nonspherical body used for calculating Jn harmonics
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_J2_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double J2, const double R_eq, const int source_index){
    const struct reb_particle source = sim->particles[source_index];
    const double G = sim->G;
    for (int i=0; i<N; i++){
        if(i == source_index){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
        const double costheta2 = dz*dz/r2;
        const double prefac = 3.*J2*R_eq*R_eq/r2/r2/r/2.;
        const double fac = 5.*costheta2-1.;

        particles[i].ax += G*source.m*prefac*fac*dx;
        particles[i].ay += G*source.m*prefac*fac*dy;
        particles[i].az += G*source.m*prefac*(fac-2.)*dz;
        particles[source_index].ax -= G*p.m*prefac*fac*dx;
        particles[source_index].ay -= G*p.m*prefac*fac*dy;
        particles[source_index].az -= G*p.m*prefac*(fac-2.)*dz;
    }
}

static void rebx_J2(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N){
    for (int i=0; i<N; i++){
        const double* const J2 = rebx_get_param_check(&particles[i], "J2", REBX_TYPE_DOUBLE);
        if (J2 != NULL){
            const double* const R_eq = rebx_get_param_check(&particles[i], "R_eq", REBX_TYPE_DOUBLE);
            if (R_eq != NULL){
                rebx_calculate_J2_force(sim, particles, N, *J2, *R_eq,i); 
            }
        }
    }
}

static void rebx_calculate_J4_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double J4, const double R_eq, const int source_index){
    const struct reb_particle source = sim->particles[source_index];
    const double G = sim->G;
    for (int i=0; i<N; i++){
        if(i == source_index){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
        const double costheta2 = dz*dz/r2;
        const double prefac = 5.*J4*R_eq*R_eq*R_eq*R_eq/r2/r2/r2/r/8.;
        const double fac = 63.*costheta2*costheta2-42.*costheta2 + 3.;

        particles[i].ax += G*source.m*prefac*fac*dx;
        particles[i].ay += G*source.m*prefac*fac*dy;
        particles[i].az += G*source.m*prefac*(fac+12.-28.*costheta2)*dz;
        particles[source_index].ax -= G*p.m*prefac*fac*dx;
        particles[source_index].ay -= G*p.m*prefac*fac*dy;
        particles[source_index].az -= G*p.m*prefac*(fac+12.-28.*costheta2)*dz;
    }
}

static void rebx_J4(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N){
    for (int i=0; i<N; i++){
        const double* const J4 = rebx_get_param_check(&particles[i], "J4", REBX_TYPE_DOUBLE);
        if (J4 != NULL){
            const double* const R_eq = rebx_get_param_check(&particles[i], "R_eq", REBX_TYPE_DOUBLE);
            if (R_eq != NULL){
                rebx_calculate_J4_force(sim, particles, N, *J4, *R_eq,i); 
            }
        }
    }
}

void rebx_gravitational_harmonics(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N){
    rebx_J2(sim, effect, particles, N);
    rebx_J4(sim, effect, particles, N);
}

static double rebx_calculate_J2_hamiltonian(struct reb_simulation* const sim, const double J2, const double R_eq, const int source_index){
    const struct reb_particle* const particles = sim->particles;
	const int _N_real = sim->N - sim->N_var;
    const struct reb_particle source = particles[source_index];
    const double G = sim->G;
    double H = 0.;
	for (int i=0;i<_N_real;i++){
		if(i == source_index){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
        const double costheta2 = dz*dz/r2;
        const double prefac = G*p.m*source.m*R_eq*R_eq/r2/r*J2;
        const double P2 = 0.5*(3.*costheta2-1.);
        H += prefac*P2;
    }		
    return H;
}

static double rebx_J2_hamiltonian(struct reb_simulation* const sim){ 
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    double Htot = 0.;
    for (int i=0; i<N_real; i++){
        const double* const J2 = rebx_get_param_check(&particles[i], "J2", REBX_TYPE_DOUBLE);
        if (J2 != NULL){
            const double* const R_eq = rebx_get_param_check(&particles[i], "R_eq", REBX_TYPE_DOUBLE);
            if (R_eq != NULL){
                Htot += rebx_calculate_J2_hamiltonian(sim, *J2, *R_eq, i); 
            }
        }
    }
    return Htot;
}

static double rebx_calculate_J4_hamiltonian(struct reb_simulation* const sim, const double J4, const double R_eq, const int source_index){
    const struct reb_particle* const particles = sim->particles;
	const int _N_real = sim->N - sim->N_var;
    const struct reb_particle source = particles[source_index];
    const double G = sim->G;
    double H = 0.;
	for (int i=0;i<_N_real;i++){
		if(i == source_index){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
        const double costheta2 = dz*dz/r2;
        const double prefac = G*p.m*source.m*R_eq*R_eq*R_eq*R_eq/r2/r2/r*J4;
        const double P4 = (35.*costheta2*costheta2 - 30.*costheta2+3.)/8.;

        H += prefac*P4;
    }		
    return H;
}

static double rebx_J4_hamiltonian(struct reb_simulation* const sim){ 
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    double Htot = 0.;
    for (int i=0; i<N_real; i++){
        const double* const J4 = rebx_get_param_check(&particles[i], "J4", REBX_TYPE_DOUBLE);
        if (J4 != NULL){
            const double* const R_eq = rebx_get_param_check(&particles[i], "R_eq", REBX_TYPE_DOUBLE);
            if (R_eq != NULL){
                Htot += rebx_calculate_J4_hamiltonian(sim, *J4, *R_eq, i); 
            }
        }
    }
    return Htot;
}

double rebx_gravitational_harmonics_hamiltonian(struct reb_simulation* const sim){ 
    double H = rebx_J2_hamiltonian(sim);
    H += rebx_J4_hamiltonian(sim);
    return H;
}
