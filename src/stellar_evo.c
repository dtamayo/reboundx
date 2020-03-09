/**
 * @file    stellar_evo.c
 * @brief   Interpolate particle parameters from external data between timesteps in the simulation.
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
 * Python Example          `StellarEvo.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/StellarEvo.ipynb>`_.
 * ======================= ===============================================
 * 
 * This interpolates particle parameter data for individual particles every timestep.
 * Set particles' ``mass_age``, ``mass_val``, and ``mass_n`` with n-sized double arrays of time-mass values to be interpolated.
 * 
 * **Effect Parameters**
 * 
 * *None*
 * 
 * **Particle Parameters**
 * 
 * Only particles with their ``mass_age``, ``mass_val``, and ``mass_n`` parameters set will have their masses affected.
 * 
 * ============================ =========== =======================================================
 * Name (C type)                Required    Description
 * ============================ =========== =======================================================
 * mass_age (double array)      Yes         Monotonic array of times in one-to-one correspondence with elements of ``mass_val``.
 * mass_val (double array)      Yes         Array of mass values in one-to-one correspondence with elements of ``mass_age``.
 * mass_n (int)                 Yes         Size of both ``mass_age`` and ``mass_val`` arrays. Mismatches will result in an invalid interpolation (``mass_n`` < actual size) or segmentation fault (``mass_n`` > actual size).
 * ============================ =========== =======================================================
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

/**
 * Given monotonic array x[0..(n-1)] and any array y[0..(n-1)], sets passed
 * particle's mass_2val[0..(n-1)] with second-order derivatives of the
 * interpolating function at the tabulated points x[i].
 * 
 * Adapted from "Numerical Recipes for C," 2nd Ed., §3.3, p. 115.
 */
void rebx_spline(struct rebx_extras* const rebx, struct reb_particle* particle, const double x[], const double y[], const int n) {
    double y2[n], u[n];
    double p, qn, sig, un;

    y2[0] = u[0] = 0.0; // lower boundary condition is set to "natural"
    for (int i=1; i<n-1; i++) {
        // the decomposition loop of the tridiagonal algorithm.
        // y2 and u are used for temporary storage of the decompsed factors.
        sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
        p = sig*y2[i-1] + 2.;
        y2[i] = (sig - 1.)/p;
        u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
        u[i] = (6.*u[i] / (x[i+1] - x[i-1]) - sig*u[i-1]) / p;
    }
    qn = un = 0.0; // upper boundary condition is set to "natural"
    y2[n-1] = (un - qn*u[n-2]) / (qn*y2[n-2] + 1.);
    for (int k=n-2; k>=0; k--) // backsubstitution loop of tridiagonal alg.
        y2[k] = y2[k]*y2[k+1] + u[k];
    rebx_set_param_pointer(rebx, particle->ap, "mass_2val", y2); // shallow copy?
}

/**
 * Given monotonic array xa[0..(n-1)], any array ya[0..(n-1)], array
 * y2a[0..(n-1)] outputted from spline above, and a value of x, this returns a
 * cubic-spline interpolated value y.
 * "Splint" comes from spl(ine)-int(erpolation).
 * 
 * Adapted from "Numerical Recipes for C," 2nd Ed., §3.3, p. 116
 */
double rebx_splint(struct rebx_extras* const rebx, struct reb_particle* particle, const double xa[], const double ya[], const double y2a[], double x) {
    int *kptr = rebx_get_param(rebx, particle->ap, "mass_klo");
    int klo;
    double h, b, a;

    if (kptr == NULL) klo = 0; // first call
    else klo = *kptr;
    // since sequential calls are in increasing order,
    // find and update place for current and future calls
    while (xa[klo+1] < x) klo++;
    rebx_set_param_int(rebx, particle->ap, "mass_klo", klo); // update klo
    h = xa[klo+1] - xa[klo];
    if (h == 0.0) {                                          // xa's must be distinct
        rebx_error(rebx, "Cubic spline run-time error...\n");
        rebx_error(rebx, "Bad xa input to routine splint\n");
        rebx_error(rebx, "...now exiting to system...\n");
        return 0;
    }
    a = (xa[klo+1] - x) / h;
    b = (x - xa[klo]) / h;
    // evaluate cubic spline
    return a*ya[klo]+b*ya[klo+1]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[klo+1])*(h*h)/6.; 
}

void rebx_stellar_evo(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt) {
    struct rebx_extras* const rebx = sim->extras;
    const int _N_real = sim->N - sim->N_var;

    for (int i=0; i<_N_real; i++) {
        struct reb_particle* const p = &sim->particles[i];
        const double* const nptr = rebx_get_param(rebx, p->ap, "mass_n");

        if (nptr != NULL) {
            const double* const xptr = rebx_get_param(rebx, p->ap, "mass_age");
            const double* const yptr = rebx_get_param(rebx, p->ap, "mass_val");

            if (xptr != NULL && yptr != NULL) {
                double* y2ptr = rebx_get_param(rebx, p->ap, "mass_2val");
                // NEED TO FIGURE OUT HOW TO DEREFERENCE
                double x[] = *xptr; // COMPILE ERROR
                double y[] = *yptr; // COMPILE ERROR
            
                if (y2ptr == NULL) {                                       // not yet splined
                    rebx_spline(rebx, p, x, y, *nptr);                     // called only once
                    y2ptr = rebx_get_param(rebx, p->ap, "mass_2val");      // immediately update
                }
                double y2[] = *y2ptr;
                double interp = rebx_splint(rebx, p, x, y, y2, sim->t+dt); // interpolate at last sim time + operator dt
                p->m = interp;
            }
            else rebx_error(rebx, "Array size set but missing values.");   // rebx_error gives meaningful err
        }
    }
    reb_move_to_com(sim);
}
