/**
 * @file    radiation_forces.c
 * @brief   Add radiation forces
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
 * $Radiation Forces$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 H. Rein, D. Tamayo
 * Implementation Paper    *In progress*
 * Based on                `Burns et al. 1979 <http://labs.adsabs.harvard.edu/adsabs/abs/1979Icar...40....1B/>`_.
 * C Example               :ref:`c_example_rad_forces_debris_disk`, :ref:`c_example_rad_forces_circumplanetary`.
 * Python Example          `Radiation_Forces_Debris_Disk.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Radiation_Forces_Debris_Disk.ipynb>`_,
 *                         `Radiation_Forces_Circumplanetary_Dust.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Radiation_Forces_Circumplanetary_Dust.ipynb>`_.
 * ======================= ===============================================
 * 
 * This applies radiation forces to particles in the simulation.  
 * It incorporates both radiation pressure and Poynting-Robertson drag.
 * Only particles whose `beta` parameter is set will feel the radiation.  
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
 * If no particles have radiation_source set, effect will assume the particle at index 0 in the particles array is the source.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * radiation_source (int)       No          Flag identifying the particle as the source of radiation.
 * beta (float)                 Yes         Ratio of radiation pressure force to gravitational force. Particles without beta set feel no radiation forces.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "reboundx.h"

static void rebx_calculate_radiation_forces(struct reb_simulation* const sim, const double c, const int source_index){
    struct reb_particle* const particles = sim->particles;
    const struct reb_particle source = particles[source_index];
    const double mu = sim->G*source.m;

    const int _N_real = sim->N - sim->N_var;
#pragma omp parallel for
    for (int i=0;i<_N_real;i++){
        
        if(i == source_index) continue;
        
        const double* beta = rebx_get_param_check(&particles[i], "beta", REBX_TYPE_DOUBLE);
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

void rebx_radiation_forces(struct reb_simulation* const sim, struct rebx_effect* const radiation_forces){ 
    double* c = rebx_get_param_check(radiation_forces, "c", REBX_TYPE_DOUBLE);
    if (c == NULL){
        reb_error(sim, "Need to set speed of light in radiation_forces effect.  See examples in documentation.\n");
    }
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    int source_found=0;
    for (int i=0; i<N_real; i++){
        if (rebx_get_param_check(&particles[i], "radiation_source", REBX_TYPE_INT) != NULL){
            source_found = 1;
            rebx_calculate_radiation_forces(sim, *c, i);
        }
    }
    if (!source_found){
        rebx_calculate_radiation_forces(sim, *c, 0);    // default source to index 0 if "radiation_source" not found on any particle
    }
}

double rebx_rad_calc_beta(const double G, const double c, const double source_mass, const double source_luminosity, const double radius, const double density, const double Q_pr){
    return 3.*source_luminosity*Q_pr/(16.*M_PI*G*source_mass*c*density*radius);   
}
double rebx_rad_calc_particle_radius(const double G, const double c, const double source_mass, const double source_luminosity, const double beta, const double density, const double Q_pr){
    return 3.*source_luminosity*Q_pr/(16.*M_PI*G*source_mass*c*density*beta);
}

