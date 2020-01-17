/**
 * @file    stark_force.c
 * @brief   Add precession forces due to tides raised on either the primary, the orbiting bodies, or both.
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
 * $Tides$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo
 * Implementation Paper    *In progress*
 * Based on                `Hut 1981 <https://ui.adsabs.harvard.edu/#abs/1981A&A....99..126H/abstract>`_.
 * C Example               :ref:`c_example_stark_force`.
 * Python Example          `TidesPrecession.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesPrecession.ipynb>`_.
 * ======================= ===============================================
 *
 * This adds precession from the tidal interactions between the particles in the simulation and the central body, both from tides raised on the primary and on the other bodies.
 * In all cases, we need to set masses for all the particles that will feel these tidal forces. After that, we can choose to include tides raised on the primary, on the "planets", or both, by setting the respective bodies' R_tides (physical radius) and k1 (apsidal motion constant, half the tidal Love number).
 * You can specify the primary with a "primary" flag.
 * If not set, the primary will default to the particle at the 0 index in the particles array.
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
 * R_tides (float)              Yes         Physical radius (required for contribution from tides raised on the body).
 * k1 (float)                   Yes         Apsidal motion constant (half the tidal Love number k2).
 * primary (int)                No          Set to 1 to specify the primary.  Defaults to treating particles[0] as primary if not set.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "reboundx.h"

// OLD VERSION W/O PARAMETER
// void rebx_stark_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
//     particles[1].ax += 0.01;
// }

void rebx_stark_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    for (int i=0; i<N; i++){
        const double* stark_acc = rebx_get_param(rebx, particles[i].ap, "stark_acc");
        if (stark_acc != NULL){
            particles[i].ax += *stark_acc;
        }
        else if (i == 1 && stark_acc == NULL) {
            printf("NULL POINTER\n");
        }
    }
}
