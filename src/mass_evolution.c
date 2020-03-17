/**
 * @file    stellar_evo.c
 * @brief   Interpolate particle parameters from passed dataset between timesteps in the simulation.
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
 * $Mass Modifications$     // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script). 
 * 
 * ======================= ===============================================
 * Authors                 Stanley A. Baronett
 * Implementation Paper    *In progress*
 * Based on                None
 * C Example               :ref:`c_example_stellar_evo`
 * Python Example          `StellarEvolution.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/StellarEvolution.ipynb>`_.
 * ======================= ===============================================
 * 
 * This interpolates particle parameter data for individual particles every timestep.
 * Set particles' ``mass_age``, ``mass_val``, and ``mass_2val`` as ``mass_n``-sized double arrays of time-mass values to be interpolated.
 * 
 * **Effect Parameters**
 * 
 * *None*
 * 
 * **Particle Parameters**
 * 
 * Only particles with their ``mass_age``, ``mass_val``, ``mass_2val``, and ``mass_n`` parameters set will have their masses affected.
 * 
 * ============================ =========== =======================================================
 * Name (C type)                Required    Description
 * ============================ =========== =======================================================
 * mass_age (double array)      Yes         Monotonic array of times in one-to-one correspondence with elements of ``mass_val``.
 * mass_val (double array)      Yes         Array of mass values in one-to-one correspondence with elements of ``mass_age``.
 * mass_2val (double array)     Yes         Uninitialized array, of size ``mass_n``, used for spline interpolation.
 * mass_n (int)                 Yes         Size of ``mass_age``, ``mass_val``, and ``mass_2val`` arrays. Mismatches will result in invalid interpolations (``mass_n`` < actual size) or segmentation faults (``mass_n`` > actual size).
 * ============================ =========== =======================================================
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "core.h"

void rebx_mass_evolution(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt) {
    struct rebx_extras* const rebx = sim->extras;

    uint32_t* hash = rebx_get_param(rebx, operator->ap, "me_hash");
    if (hash == NULL){
        return;
    }

    struct reb_particle* p = reb_get_particle_by_hash(sim, *hash);
    if (p == NULL){
        return;
    }

    const int* Nvalues = rebx_get_param(rebx, operator->ap, "me_Nvalues");
    const double* times = rebx_get_param(rebx, operator->ap, "me_times");
    const double* values = rebx_get_param(rebx, operator->ap, "me_values");

    if (Nvalues == NULL || times == NULL || values == NULL) {
        return;
    }
    
    const enum rebx_interpolation_type* interp = rebx_get_param(rebx, operator->ap, "me_interpolation");
   
    if (interp == NULL){
        rebx_set_param_int(rebx, &operator->ap, "me_interpolation", REBX_INTERPOLATION_SPLINE); // def if not set
        interp = rebx_get_param(rebx, operator->ap, "me_interpolation");
    }
    p->m = rebx_interpolate(rebx, &operator->ap, *Nvalues, sim->t, times, values, *interp);

    reb_move_to_com(sim);
}
