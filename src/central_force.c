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
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
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

static void rebx_calculate_central_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double A, const double gamma, const int source_index){
    const struct reb_particle source = particles[source_index];
    for (int i=0; i<N; i++){
        if(i == source_index){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double prefac = A*pow(r2, (gamma-1.)/2.);

        particles[i].ax += prefac*dx;
        particles[i].ay += prefac*dy;
        particles[i].az += prefac*dz;
        particles[source_index].ax -= p.m/source.m*prefac*dx;
        particles[source_index].ay -= p.m/source.m*prefac*dy;
        particles[source_index].az -= p.m/source.m*prefac*dz;
    }
}

void rebx_central_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    for (int i=0; i<N; i++){
        const double* const Acentral = rebx_get_param(sim->extras, particles[i].ap, "Acentral");
        if (Acentral != NULL){
            const double* const gammacentral = rebx_get_param(sim->extras, particles[i].ap, "gammacentral");
            if (gammacentral != NULL){
                rebx_calculate_central_force(sim, particles, N, *Acentral, *gammacentral, i); // only calculates force if a particle has both Acentral and gammacentral parameters set.
            }
        }
    }
}

static double rebx_calculate_central_force_potential(struct reb_simulation* const sim, const double A, const double gamma, const int source_index){
    const struct reb_particle* const particles = sim->particles;
	const int _N_real = sim->N - sim->N_var;
    const struct reb_particle source = particles[source_index];
    double H = 0.;
	for (int i=0;i<_N_real;i++){
		if(i == source_index){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;

        if (fabs(gamma+1.) < DBL_EPSILON){ // F propto 1/r
            H -= p.m*A*log(sqrt(r2));
        }
        else{
            H -= p.m*A*pow(r2, (gamma+1.)/2.)/(gamma+1.);
        }
    }		
    return H;
}

double rebx_central_force_potential(struct rebx_extras* const rebx){
    if (rebx->sim == NULL){
        rebx_error(rebx, ""); // rebx_error gives meaningful err
        return 0;
    }
    struct reb_simulation* sim = rebx->sim;
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    double Htot = 0.;
    for (int i=0; i<N_real; i++){
        const double* const Acentral = rebx_get_param(rebx, particles[i].ap, "Acentral");
        if (Acentral != NULL){
            const double* const gammacentral = rebx_get_param(rebx, particles[i].ap, "gammacentral");
            if (gammacentral != NULL){
                Htot += rebx_calculate_central_force_potential(sim, *Acentral, *gammacentral, i);
            }
        }
    }
    return Htot;
}

double rebx_central_force_Acentral(const struct reb_particle p, const struct reb_particle primary, const double pomegadot, const double gamma){
    struct reb_simulation* sim = p.sim;
    const double G = sim->G;
    const struct reb_orbit o = reb_tools_particle_to_orbit(G, p, primary);
    if (fabs(gamma+2.) < DBL_EPSILON){  // precession goes to 0 at r^-2, so A diverges for gamma=-2
        reb_error(sim, "Precession vanishes for force law varying as r^-2, so can't initialize Acentral from a precession rate for gamma=-2)\n");
        return 0.;
    }
    return G*primary.m*pomegadot/(1.+gamma/2.)/pow(o.d, gamma+2.)/o.n;
}
