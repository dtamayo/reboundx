/**
 * @file    gas_damping_timescale.c
 * @brief   Update orbits with prescribed timescales by directly changing orbital elements after each timestep
 * @author  Phoebe Sandhaus <pjs5535@psu.edu>
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
 * $Gas Effects$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 Phoebe Sandhaus
 * Implementation Paper    `Sandhaus et al. in prep`
 * Based on                `Dawson et al. 2016 <https://ui.adsabs.harvard.edu/abs/2016ApJ...822...54D/abstract>; Kominami & Ida 2002 <https://ui.adsabs.harvard.edu/abs/2002Icar..157...43K/abstract>`_
 * C Example               :ref:`c_example_gas_damping_timescale`
 * Python Example          `GasDampingTimescale.ipynb <https://github.com/PhoebeSandhaus/reboundx_gas_damping/tree/main/ipython_examples/GasDampingTimescale.ipynb>`_
 * ======================= ===============================================
 * 
 * This updates particles' positions and velocities between timesteps by first calculating a damping timescale for each individual particle, and then applying the timescale to damp both the eccentricity and inclination of the particle. Note: The timescale of damping should be much greater than a particle's orbital period. The damping force should also be small as compared to the gravitational forces on the particle.
 * 
 * **Effect Parameters**
 * 
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * None                         -           -
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==============================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==============================================================================
 * d_factor (double)            Yes         Depletion factor d in Equation 16 from Dawson et al. 2016; d=1 corresponds roughly to the full minimum mass solar nebula with Sigma_gas (surface gas density) = 1700 g cm^-2 at 1 AU [for d=1]; d>1 corresponds to a more depleted nebula
 * cs_coeff (double)            Yes         Sound speed coefficient; Changing the value will change assumed units; Example: If you are using the units: AU, M_sun, yr -->  cs_coeff = 0.272  # AU^(3/4) yr^-1
 * tau_coeff (double)           Yes         Timescale coefficient;   Changing the value will change assumed units; Example: If you are using the units: AU, M_sun, yr --> tau_coeff = 0.003  #  yr AU^-2
 * ============================ =========== ==============================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

static struct reb_vec3d rebx_calculate_gas_damping_timescale(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* planet, struct reb_particle* star){
    struct rebx_extras* const rebx = sim->extras;
    struct reb_orbit o = reb_orbit_from_particle(sim->G, *planet, *star);

    const double* const d_factor = rebx_get_param(rebx, planet->ap, "d_factor");
    const double* const cs_coeff = rebx_get_param(rebx, force->ap, "cs_coeff");
    const double* const tau_coeff = rebx_get_param(rebx, force->ap, "tau_coeff");
    
    struct reb_vec3d a = {0};

    if (d_factor == NULL || cs_coeff == NULL || tau_coeff == NULL){
        rebx_error(rebx, "Need to set d_factor, cs_coeff, tau_coeff parameters.  See examples in documentation.\n");
        return a;
    }

    // initialize positions and velocities
    const double dvx = planet->vx - star->vx;
    const double dvy = planet->vy - star->vy;
    const double dvz = planet->vz - star->vz;
    const double dx = planet->x-star->x;
    const double dy = planet->y-star->y;
    const double dz = planet->z-star->z;
    const double r2 = dx*dx + dy*dy + dz*dz;

    // initial semimajor axis, eccentricity, and inclination
    const double a0 = o.a;
    const double e0 = o.e;
    const double inc0 = o.inc;
    const double starMass = star->m;
    const double planetMass = planet->m;

    // eccentricity and inclination timescales from Dawson+16 Eqn 16
    double coeff;

    double vk = sqrt(sim->G*starMass/a0);
    double v = sqrt(e0*e0+inc0*inc0)*vk;
    double cs = *cs_coeff/sqrt(sqrt(a0));
    double v_over_cs = v/cs;

    if (v <= cs){
        coeff = 1.;
    }
    else {
        if (inc0 < cs/vk) {
            coeff = v_over_cs*v_over_cs*v_over_cs;
        }
        else {
            coeff = v_over_cs*v_over_cs*v_over_cs*v_over_cs;
        }
    }

    double tau_e = -*tau_coeff*(*d_factor)*a0*a0*(starMass/planetMass)*coeff;
    double tau_inc = 2.*tau_e;  // from Kominami & Ida 2002 [Eqs. 2.9 and 2.10]


    if (tau_e < INFINITY || tau_inc < INFINITY){
        const double vdotr = dx*dvx + dy*dvy + dz*dvz;
        const double prefac = 2*vdotr/r2/tau_e;
        a.x += prefac*dx;
        a.y += prefac*dy;
        a.z += prefac*dz + 2.*dvz/tau_inc;
    }
    return a;
}

void rebx_gas_damping_timescale(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    int* ptr = rebx_get_param(sim->extras, force->ap, "coordinates");
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI; // Default
    if (ptr != NULL){
        coordinates = *ptr;
    }
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_com_force(sim, force, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_gas_damping_timescale, particles, N);
}
