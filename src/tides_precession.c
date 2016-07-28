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
#include <float.h>
#include "reboundx.h"

static void rebx_calculate_tides_precession(struct reb_simulation* const sim, const int source_index){
    struct reb_particle* const particles = sim->particles;
    const struct reb_particle* source = &particles[source_index];
    const double m0 = source->m;
    const double R0 = source->r;
    double* ptr = rebx_get_param_double(source, "k2");
    double k20;
    if (ptr == NULL){
        k20 = 0.;
    }
    else{
        k20 = *ptr;
    }
    const double fac0 = k20*R0*R0*R0*R0*R0; 
    const int _N_real = sim->N - sim->N_var;
    for (int i=0;i<_N_real;i++){
        if(i == source_index) continue;
        const struct reb_particle* p = &particles[i];
        const double mratio = p->m/m0;
        if (mratio < DBL_MIN){ // m1 = 0. Continue to avoid overflow/nan
            continue;
        }
        double fac=0.; 
        fac += fac0*mratio;
        
        const double R1 = p->r;
        const double* k2 = rebx_get_param_double(p, "k2")
        if(k2 != NULL){
            fac += (*k2)*R1*R1*R1*R1*R1/mratio;
        }
        const double dx = p->x - source->x; 
        const double dy = p->y - source->y;
        const double dz = p->z - source->z;
        const double dr2 = dx*dx + dy*dy + dz*dz; 
        const double prefac = -3*sim->G*(m0 + p->m)/(dr2*dr2*dr2*dr2)*fac;

        p->ax += prefac*dx;
        p->ay += prefac*dy;
        p->az += prefac*dz;
        source->ax -= mratio*prefac*dx;
        source->ay -= mratio*prefac*dy;
        source->az -= mratio*prefac*dz;
	}
}

void rebx_tides_precession(struct reb_simulation* const sim, struct rebx_effect* const effect){ 
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    int source_found=0;
    for (int i=0; i<N_real; i++){
        if (rebx_get_param_int(&particles[i], "primary") != NULL){
            source_found = 1;
            rebx_calculate_tides_precession(sim, i);
        }
    }
    if (!source_found){
        rebx_calculate_tides_precession(sim, 0);    // default source to index 0 if "primary" not found on any particle
    }
}
