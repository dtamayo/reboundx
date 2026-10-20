/**
 * @file    galactic_tide.c
 * @brief   The effect of the galactic tidal force for a planetary system (default: the Solar System)
 * @author  Yukun Huang <yhuang.astro@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2026 Yukun Huang
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
 * Authors                 Yukun Huang
 * Implementation Paper    None
 * Based on                `Heisler & Tremaine 1986 <https://ui.adsabs.harvard.edu/abs/1986Icar...65...13H/abstract>`_.
 * C Example               :ref:`c_example_galactic_tide`
 * Python Example          `GalacticTide.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GalacticTide.ipynb>`_.
 * ======================= ===============================================
 * 
 * Adds the local Galactic tidal acceleration to bodies in a planetary system.
 * The default local stellar density, Galactic angular velocity, and unit vectors
 * are appropriate for the Solar neighborhood and are expressed in J2000 ecliptic
 * coordinates. If the simulation start time (t = 0) differs from the present
 * epoch, the user should rotate the unit vectors accordingly. For extrasolar
 * systems, the user should set the parameters to values appropriate for the
 * host star's local Galactic environment.
 * 
 * **Effect Parameters**
 *
 * All parameters are optional and are set on the force. The defaults reproduce the local
 * Galactic parameters and coordinate axes used by the implementation.
 *
 * ============================ =========== ======================================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ======================================================================================
 * gt_omega (double)            No          Galactic angular velocity in km/s/kpc (default: 26.0).
 * gt_rho (double)              No          Local density in solar masses/pc^3 (default: 0.1).
 * gt_horizontal (int)          No          Enable the horizontal tide when nonzero (default: 1).
 * gt_ex (reb_vec3d)            No          Galactic horizontal x-axis unit vector in J2000 ecliptic coordinates.
 * gt_ey (reb_vec3d)            No          Galactic horizontal y-axis unit vector in J2000 ecliptic coordinates.
 * gt_sgp (reb_vec3d)           No          South Galactic Pole unit vector in J2000 ecliptic coordinates.
 * ============================ =========== ======================================================================================
 * 
 * **Particle Parameters**
 * 
 * *None*
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

// Define 1pc in AU, according to IAU 2015 Resolution B2.
static const double GT_PC = 648000/M_PI;

// Define the default galactic coordinate axes
// Galactic anticenter (from Sun to Galactic center) direction vector in J2000 ecliptic coordinates (EX)
static const struct reb_vec3d GT_EX = {0.054873956175539, 0.9938214697787784, 0.0964766041089421};

// Sun's motion direction vector in J2000 ecliptic coordinates (EY)
static const struct reb_vec3d GT_EY = {0.4941145323285116, -0.1109899147087437, 0.8622830554839416};

// South Galactic Pole direction vector in J2000 ecliptic coordinates (SGP)
static const struct reb_vec3d GT_SGP = {8.6766334362782760e-01, 3.5360952239719772e-04, -4.9715208647342957e-01};

static double gt_get_double(
    struct rebx_extras* rebx,
    struct rebx_force* force,
    const char* name,
    const double fallback
){
    const double* value = rebx_get_param(rebx, force->ap, name);
    return value != NULL ? *value : fallback;
}

static struct reb_vec3d gt_get_vector(
    struct rebx_extras* rebx,
    struct rebx_force* force,
    const char* name,
    const struct reb_vec3d fallback
){
    const struct reb_vec3d* value =
        rebx_get_param(rebx, force->ap, name);

    return value != NULL ? *value : fallback;
}

void rebx_galactic_tide(struct reb_simulation* const sim, struct rebx_force* const tides, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;

    // If parameters are not given, use the default values
    const double OMEGA =
        gt_get_double(rebx, tides, "gt_omega", 26);
    const double RHO =
        gt_get_double(rebx, tides, "gt_rho", 0.1);
    const struct reb_vec3d EX =
        gt_get_vector(rebx, tides, "gt_ex", GT_EX);
    const struct reb_vec3d EY =
        gt_get_vector(rebx, tides, "gt_ey", GT_EY);
    const struct reb_vec3d SGP =
        gt_get_vector(rebx, tides, "gt_sgp", GT_SGP);

    // Determine if horizontal tides should be applied
    const int* horizontal_param =
        rebx_get_param(rebx, tides->ap, "gt_horizontal");

    const int horizontal =
        horizontal_param != NULL ? (*horizontal_param != 0) : 1;

    // Calculate the constant for the vertical tidal force
    const double vertical_coefficient =
        -4.0 * M_PI * sim->G * RHO / (GT_PC * GT_PC * GT_PC);

    double omega2 = 0.0;
    struct reb_vec3d ex = EX;
    struct reb_vec3d ey = EY;

    const int _N_real = sim->N - sim->N_var;

    if (horizontal){
        const double omega = OMEGA * 1.627696965e-10 * sqrt(sim->G);
        const double phase = omega * sim->t;
        const double c = cos(phase);
        const double s = sin(phase);

        omega2 = omega * omega;

        // Rotate the horizontal basis vectors according to the current phase
        ex = (struct reb_vec3d){EX.x * c + EY.x * s, EX.y * c + EY.y * s, EX.z * c + EY.z * s};
        ey = (struct reb_vec3d){-EX.x * s + EY.x * c, -EX.y * s + EY.y * c, -EX.z * s + EY.z * c};
    }

    for (int i = sim->N_active; i < _N_real; i++){
        struct reb_particle* const p = &particles[i];

        // Calculate the galactic coordinates of the particle
        const double z_gal = SGP.x * p->x + SGP.y * p->y + SGP.z * p->z;

        // Calculate the vertical tidal forces
        const double tides_z = vertical_coefficient * z_gal;

        // Update particle accelerations for the vertical tide.
        p->ax += SGP.x * tides_z;
        p->ay += SGP.y * tides_z;
        p->az += SGP.z * tides_z;

        if (horizontal){
            // Calculate the galactic coordinates of the particle in the horizontal plane
            const double x_gal = ex.x * p->x + ex.y * p->y + ex.z * p->z;
            const double y_gal = ey.x * p->x + ey.y * p->y + ey.z * p->z;

            // Calculate the horizontal tidal forces
            const double tides_x = omega2 * x_gal;
            const double tides_y = -omega2 * y_gal;

            // Update particle accelerations if horizontal tides are considered
            p->ax += ex.x * tides_x + ey.x * tides_y;
            p->ay += ex.y * tides_x + ey.y * tides_y;
            p->az += ex.z * tides_x + ey.z * tides_y;
        }
    }
}