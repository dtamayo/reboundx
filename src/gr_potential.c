/** * @file    gr_potential.c
 * @brief   Post-newtonian general relativity corrections using a simple potential that gets the pericenter precession right.
 * @author  Pengshuai (Sam) Shi, Hanno Rein, Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Pengshuai (Sam) Shi, Hanno Rein, Dan Tamayo
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
 * Authors                 H. Rein, D. Tamayo
 * Implementation Paper    *In progress*
 * Based on                `Nobili and Roxburgh 1986 <http://labs.adsabs.harvard.edu/adsabs/abs/1986IAUS..114..105N/>`_.
 * C Example               :ref:`c_example_gr`
 * Python Example          `GeneralRelativity.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GeneralRelativity.ipynb>`_.
 * ======================= ===============================================
 * 
 * This is the simplest potential you can use for general relativity.
 * It assumes that the masses are dominated by a single central body.
 * It gets the precession right, but gets the mean motion wrong by :math:`\mathcal{O}(GM/ac^2)`.  
 * It's the fastest option, and because it's not velocity-dependent, it automatically keeps WHFast symplectic.  
 * Nice if you have a single-star system, don't need to get GR exactly right, and want speed.
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
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_gr_potential(struct reb_simulation* const sim, const double C2, const int source_index){
    const int _N_real = sim->N - sim->N_var;
    const double G = sim->G;
    struct reb_particle* const particles = sim->particles;
    const struct reb_particle source = sim->particles[source_index];
    
    const double prefac1 = 6.*(G*source.m)*(G*source.m)/C2;
    for (int i=0; i<_N_real; i++){
        if(i == source_index){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double prefac = prefac1/(r2*r2);
        
        particles[i].ax -= prefac*dx;
        particles[i].ay -= prefac*dy;
        particles[i].az -= prefac*dz;
        particles[source_index].ax += p.m/source.m*prefac*dx;
        particles[source_index].ay += p.m/source.m*prefac*dy;
        particles[source_index].az += p.m/source.m*prefac*dz;
    }
}

void rebx_gr_potential(struct reb_simulation* const sim, struct rebx_effect* const gr_potential){ 
    double* c = rebx_get_param_check(gr_potential, "c", REBX_TYPE_DOUBLE);
    if (c == NULL){
        reb_error(sim, "Need to set speed of light in gr effect.  See examples in documentation.\n");
    }
    const double C2 = (*c)*(*c);
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    for (int i=0; i<N_real; i++){
        if (rebx_get_param_check(&particles[i], "gr_source", REBX_TYPE_INT) != NULL){
            rebx_calculate_gr_potential(sim, C2, i);
            return;                             // only apply effect to first gr_source found.  For multiple sources, need gr_full
        }
    }
    rebx_calculate_gr_potential(sim, C2, 0);    // gr_source not found, default to index=0
}

static double rebx_calculate_gr_potential_hamiltonian(struct reb_simulation* const sim, const double C2, const int source_index){
    const struct reb_particle* const particles = sim->particles;
	const int _N_real = sim->N - sim->N_var;
	const double G = sim->G;
    const struct reb_particle source = particles[source_index];
	const double mu = G*source.m;
    const double prefac = 3.*mu*mu/C2;
    double H = reb_tools_energy(sim);

	for (int i=0;i<_N_real;i++){
		if(i == source_index){
            continue;
        }
        struct reb_particle pi = particles[i];
        double dx = pi.x - source.x;
        double dy = pi.y - source.y;
        double dz = pi.z - source.z;
        double r2 = dx*dx + dy*dy + dz*dz;
        H -= prefac*pi.m/r2;
    }		
	
    return H;
}

double rebx_gr_potential_hamiltonian(struct reb_simulation* const sim, const struct rebx_effect* const gr_potential){
    double* c = rebx_get_param_check(gr_potential, "c", REBX_TYPE_DOUBLE);
    if (c == NULL){
        reb_error(sim, "Need to set speed of light in gr effect.  See examples in documentation.\n");
    }
    const double C2 = (*c)*(*c);
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    for (int i=0; i<N_real; i++){
        if (rebx_get_param_check(&particles[i], "gr_source", REBX_TYPE_INT) != NULL){
            return rebx_calculate_gr_potential_hamiltonian(sim, C2, i);
        }
    }
    return rebx_calculate_gr_potential_hamiltonian(sim, C2, 0);
}
