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

static void rebx_calculate_tides_precession(struct reb_simulation* const sim, const int source_index){
    struct reb_particle* const particles = sim->particles;
    struct reb_particle* const source = &particles[source_index];
    const double m0 = source->m;
    double R0 = 0.;
    double* R = rebx_get_param_check(source, "R_tides", REBX_TYPE_DOUBLE);
    if (R){
        R0 = *R;
    }
    double k10 = 0.;
    double* k1 = rebx_get_param_check(source, "k1", REBX_TYPE_DOUBLE);
    if (k1){
        k10 = *k1;
    }
    const double fac0 = k10*R0*R0*R0*R0*R0; 
    const int _N_real = sim->N - sim->N_var;
    double k1p, Rp;
    for (int i=0;i<_N_real;i++){
        if(i == source_index) continue;
        struct reb_particle* const p = &particles[i];
        const double mratio = p->m/m0;
        if (mratio < DBL_MIN){ // m1 = 0. Continue to avoid overflow/nan
            continue;
        }
        double fac=0.; 
        fac += fac0*mratio;
        
        Rp = 0.;
        R = rebx_get_param_check(p, "R_tides", REBX_TYPE_DOUBLE);
        if(R){
            Rp = *R;
        }
        k1p = 0.;
        k1 = rebx_get_param_check(p, "k1", REBX_TYPE_DOUBLE);
        if(k1){
            k1p = *k1;
        }
        
        fac += k1p*Rp*Rp*Rp*Rp*Rp/mratio;
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

void rebx_tides_precession(struct reb_simulation* const sim, struct rebx_effect* effect){
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    int source_found=0;
    for (int i=0; i<N_real; i++){
        if (rebx_get_param_check(&particles[i], "primary", REBX_TYPE_INT) != NULL){
            source_found = 1;
            rebx_calculate_tides_precession(sim, i);
        }
    }
    if (!source_found){
        rebx_calculate_tides_precession(sim, 0);    // default source to index 0 if "primary" not found on any particle
    }
}

static double rebx_calculate_tides_precession_hamiltonian(struct reb_simulation* const sim, const int source_index){
    struct reb_particle* const particles = sim->particles;
    struct reb_particle* const source = &particles[source_index];
    const double m0 = source->m;
    double R0 = 0.;
    double* R = rebx_get_param_check(source, "R_tides", REBX_TYPE_DOUBLE);
    if (R){
        R0 = *R;
    }
    double k10 = 0.;
    double* k1 = rebx_get_param_check(source, "k1", REBX_TYPE_DOUBLE);
    if (k1){
        k10 = *k1;
    }
    const double fac0 = k10*R0*R0*R0*R0*R0; 
    const int _N_real = sim->N - sim->N_var;
    double H=0.;
    double Rp, k1p;
    for (int i=0;i<_N_real;i++){
        if(i == source_index) continue;
        struct reb_particle* const p = &particles[i];
        const double mratio = p->m/m0;
        if (mratio < DBL_MIN){ // m1 = 0. Continue to avoid overflow/nan
            continue;
        }
        double fac=0.; 
        fac += fac0*mratio;
        
        Rp = 0.;
        R = rebx_get_param_check(p, "R_tides", REBX_TYPE_DOUBLE);
        if(R){
            Rp = *R;
        }
        k1p = 0.;
        k1 = rebx_get_param_check(p, "k1", REBX_TYPE_DOUBLE);
        if(k1){
            k1p = *k1;
        }
        
        fac += k1p*Rp*Rp*Rp*Rp*Rp/mratio;
        
        const double dx = p->x - source->x; 
        const double dy = p->y - source->y;
        const double dz = p->z - source->z;
        const double dr2 = dx*dx + dy*dy + dz*dz; 
        H += -3./6.*sim->G*(m0 + p->m)*p->m/(dr2*dr2*dr2)*fac;
	}
    return H;
}

double rebx_tides_precession_hamiltonian(struct reb_simulation* const sim){
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    int source_found=0;
    double H=0.;
    for (int i=0; i<N_real; i++){
        if (rebx_get_param_check(&particles[i], "primary", REBX_TYPE_INT) != NULL){
            source_found = 1;
            H = rebx_calculate_tides_precession_hamiltonian(sim, i);
        }
    }
    if (!source_found){
        H = rebx_calculate_tides_precession_hamiltonian(sim, 0);    // default source to index 0 if "primary" not found on any particle
    }
    return H;
}
