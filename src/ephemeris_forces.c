/** * @file central_force.c
 * @brief   A general central force.
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
 * $Central Force$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo
 * Implementation Paper    *In progress*
 * Based on                None
 * C Example               :ref:`c_example_central_force`
 * Python Example          `CentralForce.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/CentralForce.ipynb>`_.
 * ======================= ===============================================
 * 
 * Adds a general central acceleration of the form a=Acentral*r^gammacentral, outward along the direction from a central particle to the body.
 * Effect is turned on by adding Acentral and gammacentral parameters to a particle, which will act as the central body for the effect,
 * and will act on all other particles.
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
 * Acentral (double)             Yes         Normalization for central acceleration.
 * gammacentral (double)         Yes         Power index for central acceleration.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

static void ephem(const int i, const double t, double* const m, double* const x, double* const y, double* const z){
    const double n = 1.;
    const double mu = 1.e-3;
    const double m0 = 1.-mu;
    const double m1 = mu; 
    if (i==0){
        *m = m0;
        const double mfac = -m1/(m0+m1);
        *x = mfac*cos(n*t);
        *y = mfac*sin(n*t);
        *z = 0.;
    }

    if (i==1){
        *m = m1;
        const double mfac = m0/(m0+m1);
        *x = mfac*cos(n*t);
        *y = mfac*sin(n*t);
        *z = 0.;
    }
}

void rebx_ephemeris_forces(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    const int* const N_ephem = rebx_get_param(sim->extras, force->ap, "N_ephem");
    if (N_ephem == NULL){
        fprintf(stderr, "REBOUNDx Error: Need to set N_ephem for ephemeris_forces\n");
        return;
    }

    const double G = sim->G;
    const double t = sim->t;
    double m, x, y, z;
    for (int i=0; i<*N_ephem; i++){
        ephem(i, t, &m, &x, &y, &z);
        for (int j=0; j<N; j++){
            const double dx = particles[j].x - x;
            const double dy = particles[j].y - y;
            const double dz = particles[j].z - z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz);
            const double prefac = G*m/(_r*_r*_r);
            //fprintf(stderr, "%e, %e, %e, %e\n", m, x, y, z);
            particles[j].ax -= prefac*dx;
            particles[j].ay -= prefac*dy;
            particles[j].az -= prefac*dz;
        }
    }
}
