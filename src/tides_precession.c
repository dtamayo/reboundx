/**
 * @file    tides_precession.c
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
 * C Example               :ref:`c_example_tides_precession`.
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

static void rebx_calculate_tides(struct reb_particle* source, struct reb_particle* target, const double G, const double k1, const double tau){
    const double ms = source->m;
    const double mt = target->m;
    const double Rt = target->r;

    const double mratio = ms/mt; // have already checked for 0 and inf
    const double fac = mratio*k1*Rt*Rt*Rt*Rt*Rt; 
    
    const double dx = target->x - source->x; 
    const double dy = target->y - source->y;
    const double dz = target->z - source->z;
    const double dr2 = dx*dx + dy*dy + dz*dz; 
    const double prefac = -3*G/(dr2*dr2*dr2*dr2)*fac;

    target->ax += prefac*ms*dx;
    target->ay += prefac*ms*dy;
    target->az += prefac*ms*dz;
    source->ax -= prefac*mt*dx;
    source->ay -= prefac*mt*dy;
    source->az -= prefac*mt*dz;
}


void rebx_tides_precession(struct reb_simulation* const sim, struct rebx_force* const tides_prec, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    const double G = sim->G;

    // Calculate tides raised on star
    struct reb_particle* target = &particles[0];// assumes nearly Keplerian motion around a single primary (particles[0])
    if (target->m == 0){                        // nothing makes sense if primary has no mass
        return;
    }
    double* k1 = rebx_get_param(rebx, target->ap, "tctl_k1");
    if (k1 != NULL && target->r != 0){  // tides on star only nonzero if k1 and finite size are set
        // We don't require time lag tau to be set. Might just want conservative piece of tidal potential
        double tau = 0.;
        double* tauptr = rebx_get_param(rebx, target->ap, "tctl_tau");
        if (tauptr){
            tau = *tauptr;
        }
        for (int i=1; i<N; i++){
            struct reb_particle* source = &particles[i]; // planet raising the tides on the star
            if (source->m == 0){
                continue;
            }
            rebx_calculate_tides(source, target, G, *k1, tau);
        }
    }

    // Calculate tides raised on the planets
    struct reb_particle* source = &particles[0]; // Source is always the star (no planet-planet tides)
    for (int i=1; i<N; i++){
        struct reb_particle* target = &particles[i]; 
        double* k1 = rebx_get_param(rebx, target->ap, "tctl_k1");
        if (k1 == NULL || target->r == 0 || target->m == 0){
            continue;
        }
        double tau = 0.;
        double* tauptr = rebx_get_param(rebx, target->ap, "tctl_tau");
        if (tauptr){
            tau = *tauptr;
        }
        rebx_calculate_tides(source, target, G, *k1, tau);
    }
}

static double rebx_calculate_tides_potential(struct reb_particle* source, struct reb_particle* target, const double G, const double k1){
    const double ms = source->m;
    const double mt = target->m;
    const double Rt = target->r;

    const double mratio = ms/mt; // have already checked for 0 and inf
    const double fac = mratio*k1*Rt*Rt*Rt*Rt*Rt; 
    
    const double dx = target->x - source->x; 
    const double dy = target->y - source->y;
    const double dz = target->z - source->z;
    const double dr2 = dx*dx + dy*dy + dz*dz; 
    
    return -1./2.*G*ms*mt/(dr2*dr2*dr2)*fac;
}

double rebx_tides_precession_potential(struct rebx_extras* const rebx){
    if (rebx->sim == NULL){
        rebx_error(rebx, ""); // rebx_error gives meaningful err
        return 0;
    }
    struct reb_simulation* const sim = rebx->sim;
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    const double G = sim->G;
    double H=0.;

    // Calculate tides raised on star
    struct reb_particle* target = &particles[0];// assumes nearly Keplerian motion around a single primary (particles[0])
    if (target->m == 0){                        // No potential with massless primary
        return 0.;
    }
    double* k1 = rebx_get_param(rebx, target->ap, "tctl_k1");
    if (k1 != NULL && target->r != 0){  // tides on star only nonzero if k1 and finite size are set
        for (int i=1; i<N_real; i++){
            struct reb_particle* source = &particles[i]; // planet raising the tides on the star
            if (source->m == 0){
                continue;
            }
            H += rebx_calculate_tides_potential(source, target, G, *k1);
        }
    }

    // Calculate tides raised on the planets
    struct reb_particle* source = &particles[0]; // Source is always the star (no planet-planet tides)
    for (int i=1; i<N_real; i++){
        struct reb_particle* target = &particles[i]; 
        double* k1 = rebx_get_param(rebx, target->ap, "tctl_k1");
        if (k1 == NULL || target->r == 0 || target->m == 0){
            continue;
        }
        H += rebx_calculate_tides_potential(source, target, G, *k1);
    }

    return H;
}
