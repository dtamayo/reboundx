/**
 * @file    tides_drag.c
 * @brief   Add drag forces due to slowly rotating tides raised on the primary body.
 * @author  Stanley A. Baronett <stanley.a.baronett@gmail.com>
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
 * Authors                 S.A. Baronett
 * Implementation Paper    *In progress*
 * Based on                `Schroder & Smith 2008 <https://arxiv.org/abs/0801.4031>`_.
 * C Example               :ref:`c_example_tides_drag`.
 * Python Example          `TidesDrag.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesDrag.ipynb>`_.
 * ======================= ===============================================
 *
 * This adds drag from the tidal interactions between the particles in the simulation and the central body, from slowly rotating tides raised on the primary body.
 * The effect assumes all affected particles follow co-planar, circular orbits with respect to the primary's rotational pole.
 * In all cases, we need to set masses for all the particles that will feel these tidal forces. After that, we can choose to include tides raised on the primary by setting its R_tides (physical radius), luminosity (in units consistent with G), and tides_lambda2 (a coefficient that depends on the properties of the primary's convective envelope).
 * You can specify the primary with a "primary" flag and its rotational angular velocity with the "Omega" parameter.
 * If not set, the primary will default to the particle at the 0 index in the particles array and its rotational angular velocity will default to 0.
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
 * R_tides (float)              Yes         Primary physical radius (required for contribution from tides raised on the body).
 * luminosity (float)           Yes         Primary luminosity (in units consistent with G).
 * tides_lambda2 (float)        Yes         Coefficient that depends on the properties of the primary convective envelope.
 * primary (int)                No          Set to 1 to specify the primary.  Defaults to treating particles[0] as primary if not set.
 * Omega (float)                No          Primary rotational angular velocity of the primary.  Defaults to 0 if not set.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "reboundx.h"

static void rebx_calculate_tides_drag(struct rebx_extras* const rebx, struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const int source_index){
    struct reb_particle* const source = &particles[source_index];
    double* R  = rebx_get_param(rebx, source->ap, "R_tides");
    double* L  = rebx_get_param(rebx, source->ap, "luminosity");
    double* l2 = rebx_get_param(rebx, source->ap, "tides_lambda2");
    double* O  = rebx_get_param(rebx, source->ap, "tides_Omega");
    
    if      (R == NULL)  rebx_error(rebx, "REBOUNDx Error: Primary \"R_tides\" (physical radius) required. Ignoring effect.");
    else if (L == NULL)  rebx_error(rebx, "REBOUNDx Error: Primary \"luminosity\" required. Ignoring effect.");   
    else if (l2 == NULL) rebx_error(rebx, "REBOUNDx Error: Primary \"lambda2\" required. Ignoring effect.");
    else{
        struct reb_orbit po;
        double m, dr, vmag, omega, q, rratio;
        double R0 = *R;                 
        double L0 = *L;                 
        double M0 = source->m;          
        double t_f = cbrt(M0*R0*R0/L0); // convective friction time (Zahn 1989, Eq.15)
        double lambda2 = *l2;           
        double Omega = 0;               
        if (O != NULL) Omega = *O;
        
        for (int i=0; i<N; i++){
            if (i == source_index) continue;
            po = reb_tools_particle_to_orbit(sim->G, particles[i], particles[0]);
            dr = po.d;
            vmag = po.v;
            m = particles[i].m;
            q = m/M0;                   // mass ratio
            rratio = R0/dr;             // ratio of primary physical radius to orbital radius
            omega = po.n;

            if (q < DBL_MIN) continue;  // m1 = 0. Continue to avoid overflow/nan

            // Equation (4) of Schroder & Smith (2008):
            const double torque = 6.*(lambda2/t_f)*q*q*M0*R0*R0*rratio*rratio*rratio*rratio*rratio*rratio*(Omega - omega);

            const double prefac = torque/m/dr/vmag;

            particles[i].ax += prefac/(1-q)*(particles[i].vx - source->vx);
            particles[i].ay += prefac/(1-q)*(particles[i].vy - source->vy);
            particles[i].az += prefac/(1-q)*(particles[i].vz - source->vz);
            source->ax -= prefac*q/(1-q)*(particles[i].vx - source->vx);
            source->ay -= prefac*q/(1-q)*(particles[i].vy - source->vy);
            source->az -= prefac*q/(1-q)*(particles[i].vz - source->vz);
        }
    }
}

void rebx_tides_drag(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    _Bool source_found = 0;
    
    for (int i=0; i<N; i++){
        if (rebx_get_param(rebx, particles[i].ap, "tides_primary") != NULL){
            source_found = 1;
            rebx_calculate_tides_drag(rebx, sim, particles, N, i);
        }
    }
    if (!source_found) rebx_calculate_tides_drag(rebx, sim, particles, N, 0); // default source to index 0 if "tides_primary" not found on any particle
}

/*
CONVENIENCE FUNCTIONS
*/
double rebx_calculate_tides_drag_lambda2(const double alpha){
    return cbrt(0.019*alpha*alpha*alpha*alpha);
}

double rebx_calculate_tides_drag_t_f(const double primary_mass, const double primary_radius, const double source_luminosity){
    return cbrt(primary_mass*primary_radius*primary_radius/source_luminosity);
}

double rebx_calculate_tides_drag_q(const double secondary_mass, const double primary_mass){
    return secondary_mass/primary_mass;
}

double rebx_calculate_tides_drag_Gamma(const double lambda2, const double t_f, const double q, const double primary_mass, const double primary_radius, const double orbital_radius, const double Omega, const double omega){
    const double rratio = primary_radius / orbital_radius;
    return 6.*(lambda2/t_f)*q*q*primary_mass*primary_radius*primary_radius*rratio*rratio*rratio*rratio*rratio*rratio*(Omega - omega);
}
