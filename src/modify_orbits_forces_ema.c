/**
 * @file    modify_orbits_forces_ema.c
 * @brief   Update orbital elements with prescribed timescales using forces.
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
 * $Orbit Modifications$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo, H. Rein
 * Implementation Paper    `Kostov et al., 2016 <https://ui.adsabs.harvard.edu/abs/2016ApJ...832..183K/abstract>`_.
 * Based on                `Papaloizou & Larwood 2000 <http://labs.adsabs.harvard.edu/adsabs/abs/2000MNRAS.315..823P/>`_.
 * C Example               :ref:`c_example_modify_orbits`
 * Python Example          `Migration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Migration.ipynb>`_
 *                         `EccAndIncDamping.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/EccAndIncDamping.ipynb>`_.
 * ======================= ===============================================
 * 
 * This applies physical forces that orbit-average to give exponential growth/decay of the semimajor axis, eccentricity and inclination.
 * The eccentricity damping keeps the angular momentum constant (corresponding to `p=1` in modify_orbits_direct), which means that eccentricity damping will induce some semimajor axis evolution.
 * Additionally, eccentricity/inclination damping will induce pericenter/nodal precession.
 * Both these effects are physical, and the method is more robust for strongly perturbed systems.
 * 
 * **Effect Parameters**
 *
 * If coordinates not, defaults to using Jacobi coordinates.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * coordinates (enum)           No          Type of elements to use for modification (Jacobi, barycentric or particle).
 *                                          See the examples for usage.
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * One can pick and choose which particles have which parameters set.  
 * For each particle, any unset parameter is ignored.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * tau_a (double)               No          Semimajor axis exponential growth/damping timescale
 * tau_e (double)               No          Eccentricity exponential growth/damping timescale
 * tau_inc (double)             No          Inclination axis exponential growth/damping timescale
 * p (double)                   No          Amount of angular momentun preserved
 * tdisk (double)               No          Characteristic lifetime of the disk
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

static struct reb_vec3d rebx_calculate_modify_orbits_forces_ema(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p, struct reb_particle* source){
    double invtau_a = 0.0;
    double tau_e = INFINITY;
    double tau_inc = INFINITY;
    
    const double* const tau_a_ptr = rebx_get_param(sim->extras, p->ap, "tau_a");
    const double* const tau_e_ptr = rebx_get_param(sim->extras, p->ap, "tau_e");
    const double* const tau_inc_ptr = rebx_get_param(sim->extras, p->ap, "tau_inc");

    double pval = 1.;
    double tdisk = INFINITY;

    const double* const p_ptr = rebx_get_param(sim->extras, force->ap, "p");
    const double* const tdisk_ptr = rebx_get_param(sim->extras, force->ap, "tdisk");

    const double dvx = p->vx - source->vx;
    const double dvy = p->vy - source->vy;
    const double dvz = p->vz - source->vz;
    const double dx = p->x - source->x;
    const double dy = p->y - source->y;
    const double dz = p->z - source->z;
    const double r2 = dx*dx + dy*dy + dz*dz;
    
    if(tau_a_ptr != NULL){
        invtau_a = 1.0/(*tau_a_ptr);
    }
    if(tau_e_ptr != NULL){
        tau_e = *tau_e_ptr;
    }
    if(tau_inc_ptr != NULL){
        tau_inc = *tau_inc_ptr;
    }
    if(p_ptr != NULL){
        pval = *p_ptr;
    }
    if(tdisk_ptr != NULL){
        tdisk = *tdisk_ptr;
    }
    
    if(pval > 0. && tau_e < INFINITY){
        int err = 0;
        struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *p, *source, &err);
        const double e1 = o.e;
        invtau_a += 2 * pval * e1 * e1 / tau_e;
    }
    
    struct reb_vec3d a = {0};
    a.x = dvx * invtau_a / 2.;
    a.y = dvy * invtau_a / 2.;
    a.z = dvz * invtau_a / 2.;

    if(tau_e < INFINITY){
        const double vdotr = dx*dvx + dy*dvy + dz*dvz;
        const double prefac = 2 * vdotr / r2 / tau_e;
        a.x += prefac * dx;
        a.y += prefac * dy;
        a.z += prefac * dz;
    }
    
    if(tau_inc < INFINITY){
        a.z += 2 * dvz / tau_inc;
    }
    
    if(tdisk < INFINITY){
        double fac_stokes = 0.5 * (1. + tanh(10 * (1. - sim->t / tdisk)));
        a.x *= fac_stokes;
        a.y *= fac_stokes;
        a.z *= fac_stokes;
    }
    
    return a;
}

void rebx_modify_orbits_forces_ema(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    int* ptr = rebx_get_param(sim->extras, force->ap, "coordinates");
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI; // Default
    if (ptr != NULL){
        coordinates = *ptr;
    }
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_com_force(sim, force, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_modify_orbits_forces_ema, particles, N);
}
