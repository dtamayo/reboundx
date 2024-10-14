/** * @file j2.c
 * @brief   Adding J2 to Phoebe.
 * @author  Miroslav Broz <miroslav.broz@email.cz>
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
 * Authors                 M. Broz
 * Implementation Paper    `Fabrycky, 2010 <https://ui.adsabs.harvard.edu/abs/2010exop.book..217F/abstract>`_. 
 * Based on                None
 * ======================= ===============================================
 * 
 * Adds azimuthally symmetric gravitational harmonics (J2 = -C20), or oblateness. Alternative implementation, assuming that each pair of bodies has spin poles perpendicular to its orbital plane. This approximation is suitable for stellar systems.
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

/*
 * From Fabrycky (2010), Eq. (4):
 *
 * f_R = -1/2 k_L Omega_rot^2 R^5 / r^5 \vec r
 *
 * J2 = -C20
 * n = sqrt(GM/R^3)
 * k_L = -C20 (n/Omega_rot)^2
 *
 * Note: This includes only the radial component of acceleration.
 *
 */

void rebx_j2_phoebe(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){

    const double G = sim->G;
    struct rebx_extras* const rebx = sim->extras;

    for (int i=0; i<N; i++){
        const double* const J2 = rebx_get_param(rebx, particles[i].ap, "J2");
        if (J2 == NULL){
            continue;
        }
        const double* const R_eq = rebx_get_param(rebx, particles[i].ap, "R_eq");
        if (R_eq == NULL){
            continue;
        }
        const struct reb_particle pi = particles[i];

        for (int j=0; j<N; j++){
            if (j == i){
                continue;
            }
            const struct reb_particle pj = particles[j];

            const double dx = pj.x - pi.x;
            const double dy = pj.y - pi.y;
            const double dz = pj.z - pi.z;
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double r = sqrt(r2);
            double fac = -0.5*(*J2)*G*pi.m*(*R_eq)*(*R_eq)/r2/r2/r;

            particles[j].ax += fac*dx;
            particles[j].ay += fac*dy;
            particles[j].az += fac*dz;
         
            fac *= pj.m/pi.m;
         
            particles[i].ax -= fac*dx;
            particles[i].ay -= fac*dy;
            particles[i].az -= fac*dz;
        }
    }
}


