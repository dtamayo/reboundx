/**
 * @file    tides_synchronous_ecc_damping.c
 * @brief   Adds the tidal eccentricity damping that does not depend on spins and exists even for synchronous bodies
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
 * C Example               :ref:``.
 * Python Example          ` <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesPrecession.ipynb>`_.
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

#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_tides_synchronous_ecc_damping(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const int source_index){
    const struct reb_particle source = particles[source_index];
    for(int i=0; i<N; i++){
        if(i==source_index){
            continue;
        }
        const double* const tau_e_ptr = rebx_get_param_check(&particles[i], "tides_synchronous_tau_e", REBX_TYPE_DOUBLE);
        if(tau_e_ptr != NULL){
            const double tau_e = *tau_e_ptr;
            const struct reb_particle p = particles[i];
            const double dvx = p.vx - source.vx;
            const double dvy = p.vy - source.vy;
            const double dvz = p.vz - source.vz;

            const double dx = p.x-source.x;
            const double dy = p.y-source.y;
            const double dz = p.z-source.z;
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double vdotr = dx*dvx + dy*dvy + dz*dvz;
            const double prefac = 2.*vdotr/r2/tau_e;
            const double prefacp = prefac*source.m/(p.m+source.m);
            const double prefacprimary = prefac*p.m/(p.m+source.m);

            particles[i].ax += prefacp*dx;
            particles[i].ay += prefacp*dy;
            particles[i].az += prefacp*dz;
            
            particles[0].ax -= prefacprimary*dx;
            particles[0].ay -= prefacprimary*dy;
            particles[0].az -= prefacprimary*dz;
        }
    }
}

void rebx_tides_synchronous_ecc_damping(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N){
    int source_found=0;
    for (int i=0; i<N; i++){
        if (rebx_get_param_check(&particles[i], "tides_primary", REBX_TYPE_INT) != NULL){
            source_found = 1;
            rebx_calculate_tides_synchronous_ecc_damping(sim, particles, N, i);
        }
    }
    if (!source_found){
        rebx_calculate_tides_synchronous_ecc_damping(sim, particles, N, 0);    // default source to index 0 if "tides_primary" not found on any particle
    }
}
