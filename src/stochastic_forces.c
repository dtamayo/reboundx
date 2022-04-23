/**
 * @file    stochastic_forces.c
 * @brief   Add stochastic forces
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE
 * Copyright (c) 2022 Hanno Rein, Dan Tamayo
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
 * $Stochastic Forces$     
 *
 * ======================= ===============================================
 * Authors                 H. Rein
 * Based on                `Rein and Papaloizou 2009 <https://ui.adsabs.harvard.edu/abs/2009A%26A...497..595R/abstract>`_.
 * C Example               :ref:`c_example_stochastic_forces`, :ref:`c_example_stochastic_forces`.
 * Python Example          `Stochastic_Forces.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Stochastic_Forces.ipynb>`_,
 * ======================= ===============================================
 * 
 * This applies stochastic forces to particles in the simulation.  
 * 
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 *
 * All particles which have the field stochastic_D set, will experience stochastic forces.
 * The particle with index 0 cannot experience stochastic forces.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * stochastic_D (float)         Yes         Diffusion coefficient
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "reboundx.h"

static void rebx_calculate_stochastic_forces(struct rebx_extras* const rebx, struct reb_simulation* const sim, const double c, const int source_index, struct reb_particle* const particles, const int N){
    const struct reb_particle source = particles[source_index];
    const double mu = sim->G*source.m;

    for (int i=0;i<N;i++){
        
        if(i == source_index) continue;
        
        const double* beta = rebx_get_param(rebx, particles[i].ap, "beta");
        if(beta == NULL) continue; // only particles with beta set feel radiation forces
        
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x; 
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double dr = sqrt(dx*dx + dy*dy + dz*dz); // distance to star
        
        const double dvx = p.vx - source.vx;
        const double dvy = p.vy - source.vy;
        const double dvz = p.vz - source.vz;
        const double rdot = (dx*dvx + dy*dvy + dz*dvz)/dr; // radial velocity
        const double a_rad = *beta*mu/(dr*dr);

        // Equation (5) of Burns, Lamy & Soter (1979)

        particles[i].ax += a_rad*((1.-rdot/c)*dx/dr - dvx/c);
        particles[i].ay += a_rad*((1.-rdot/c)*dy/dr - dvy/c);
        particles[i].az += a_rad*((1.-rdot/c)*dz/dr - dvz/c);
	}
}

void rebx_stochastic_forces(struct reb_simulation* const sim, struct rebx_force* const radiation_forces, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    double* c = rebx_get_param(rebx, radiation_forces->ap, "c");
    if (c == NULL){
        reb_error(sim, "Need to set speed of light in radiation_forces effect.  See examples in documentation.\n");
    }
    
    int source_found=0;
    for (int i=0; i<N; i++){
        if (rebx_get_param(rebx, particles[i].ap, "radiation_source") != NULL){
            source_found = 1;
            rebx_calculate_radiation_forces(rebx, sim, *c, i, particles, N);
        }
    }
    if (!source_found){
        rebx_calculate_radiation_forces(rebx, sim, *c, 0, particles, N);    // default source to index 0 if "radiation_source" not found on any particle
    }
}

