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
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_. 
 * Based on                None
 * C Example               :ref:`c_example_J2`
 * Python Example          `J2.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/J2.ipynb>`_.
 * ======================= ===============================================
 * 
 * Adds azimuthally symmetric gravitational harmonics (J2, J4) to bodies in the simulation. Current implementation assumes everything is planar, i.e. spin pole of body aligned with z axis of simulation.
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
 * R_eq (double)                No          Equatorial radius of nonspherical body used for calculating Jn harmonics
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_lense_thirring(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){

    double Omega_x, Omega_y, Omega_z, v_x, v_y, v_z;

    Omega_x=0;
    Omega_y=0;
    Omega_z=0.01;
    
    v_x = particles[1].vx;
    v_y = particles[1].vy;
    v_z = particles[1].vz;

   particles[1].ax += 2.*Omega_y*v_z - Omega_z*v_y;
   particles[1].ay += 2.*Omega_z*v_x - Omega_x*v_z;
   particles[1].az += 2.*Omega_x*v_y - Omega_y*v_x;

}

/*
void rebx_gravitational_harmonics(struct reb_simulation* const sim, struct rebx_force* const gh, struct reb_particle* const particles, const int N){

    rebx_J2(sim->extras, sim, gh, particles, N);
}
*/

/*
static void rebx_calculate_J2_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double J2, const double R_eq, const int source_index){
    const struct reb_particle source = particles[source_index];
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

static void rebx_J2(struct rebx_extras* const rebx, struct reb_simulation* const sim, struct rebx_force* const gh, struct reb_particle* const particles, const int N){
    for (int i=0; i<N; i++){
        const double* const J2 = rebx_get_param(rebx, particles[i].ap, "J2");
        if (J2 != NULL){
            const double* const R_eq = rebx_get_param(rebx, particles[i].ap, "R_eq");
            if (R_eq != NULL){
                rebx_calculate_J2_force(sim, particles, N, *J2, *R_eq,i); 
            }
        }
    }
}

*/
