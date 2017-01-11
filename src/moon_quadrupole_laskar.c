/** * @file moon_quadrupole_laskar.c
 * @brief   Models Earth's moon as a quadrupole around the Earth interacting with the Sun
 * @author  Aleksandar Rachkov, Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Aleksandar Rachkov, Dan Tamayo, Hanno Rein
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
 * $Multipoles$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 A. Rachkov, D. Tamayo
 * Implementation Paper    *In progress*
 * Based on                None
 * C Example               :ref:`c_example_moon_quadrupole_laskar`
 * Python Example          `MoonQuadrupoleLaskar.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/MoonQuadrupoleLaskar.ipynb>`_.
 * ======================= ===============================================
 * 
 * Add description here 
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
 * Aradial (double)             Yes         Normalization for central acceleration.
 * gammaradial (double)         Yes         Power index for central acceleration.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_force(struct reb_simulation* const sim, const double A, const double gamma, const int source_index){
    const int _N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    const struct reb_particle source = sim->particles[0];
    const struct reb_particle p = particles[i];
    const double dx = p.x - source.x;
    const double dy = p.y - source.y;
    const double dz = p.z - source.z;
    const double r2 = dx*dx + dy*dy + dz*dz;

    // Calculate A. Put in right value of gamma
    const double prefac = A*pow(r2, (gamma-1.)/2.);

    particles[i].ax += prefac*dx;
    particles[i].ay += prefac*dy;
    particles[i].az += prefac*dz;
    particles[0].ax -= p.m/source.m*prefac*dx;
    particles[0].ay -= p.m/source.m*prefac*dy;
    particles[0].az -= p.m/source.m*prefac*dz;
}

void rebx_moon_quadrupole_laskar(struct reb_simulation* const sim, struct rebx_effect* const effect){ 
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    for (int i=1; i<N_real; i++){
        const double* const m_moon_mql = rebx_get_param_check(&particles[i], "m_moon_mql", REBX_TYPE_DOUBLE);
        if (m_moon_mql != NULL){
            const double* const a0_mql = rebx_get_param_check(&particles[i], "a0_mql", REBX_TYPE_DOUBLE);
            if (a0_mql != NULL){
                // keep chaining ifs to get all the parameters
                rebx_calculate_force(sim, *m_moon_mql, *a0, *a1, *a2, *alpha, i); // only calculates force if all parameters set
            }
        }
    }
}

static double rebx_calculate_hamiltonian(struct reb_simulation* const sim, const double A, const double gamma, const int source_index){
    const struct reb_particle* const particles = sim->particles;
	const int _N_real = sim->N - sim->N_var;
    const struct reb_particle source = particles[0];
    double H = 0.;
    const struct reb_particle p = particles[i];
    const double dx = p.x - source.x;
    const double dy = p.y - source.y;
    const double dz = p.z - source.z;
    const double r2 = dx*dx + dy*dy + dz*dz;

    H -= p.m*A*pow(r2, (gamma+1.)/2.)/(gamma+1.); // figure this out
   
    return H;
}

double rebx_moon_quadrupole_laskar_hamiltonian(struct reb_simulation* const sim){ 
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    double Htot = 0.;
    for (int i=1; i<N_real; i++){
        // copy paste the chain of if statements from above
        const double* const Acentral = rebx_get_param_check(&particles[i], "Acentral", REBX_TYPE_DOUBLE);
        if (Acentral != NULL){
            const double* const gammacentral = rebx_get_param_check(&particles[i], "gammacentral", REBX_TYPE_DOUBLE);
            if (gammacentral != NULL){
                Htot += rebx_calculate_hamiltonian(sim, *m_moon_mql, *a0, *a1, *a2, *alpha, i);
            }
        }
    }
    return Htot;
}
