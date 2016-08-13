/**
 * @file    modify_orbits_direct.c
 * @brief   Update orbital with prescribed timescales by directly changing orbital elements after each timestep.
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
 * Authors                 D. Tamayo
 * Implementation Paper    *In progress*
 * Based on                `Lee & Peale 2002 <http://labs.adsabs.harvard.edu/adsabs/abs/2002ApJ...567..596L/>`_. 
 * C Example               :ref:`c_example_modify_orbits`
 * Python Example          `Migration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Migration.ipynb>`_,
 *                         `EccAndIncDamping.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/EccAndIncDamping.ipynb>`_.
 * ======================= ===============================================
 * 
 * This updates particles' positions and velocities between timesteps to achieve the desired changes to the osculating orbital elements (exponential growth/decay for a, e, inc, linear progression/regression for Omega/omega.
 * This nicely isolates changes to particular osculating elements, making it easier to interpret the resulting dynamics.  
 * One can also adjust the coupling parameter `p` between eccentricity and semimajor axis evolution, as well as whether the damping is done on Jacobi, barycentric or heliocentric elements.
 * Since this method changes osculating (i.e., two-body) elements, it can give unphysical results in highly perturbed systems.
 * 
 * **Effect Parameters**
 *
 * If p is not set, it defaults to 1.  If coordinates not set, defaults to using Jacobi coordinates.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * p (double)                   No          Coupling parameter between eccentricity and semimajor axis evolution
 *                                          (see Deck & Batygin 2015). `p=0` corresponds to no coupling, `p=1` to
 *                                          eccentricity evolution at constant angular momentum.
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
 * tau_Omega (double)           No          Period of linear nodal precession/regression
 * tau_omega (double)           No          Period of linear apsidal precession/regression
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

static struct reb_particle rebx_calculate_modify_orbits_direct(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* p, struct reb_particle* primary){
    int err=0;
    struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *p, *primary, &err);
    if(err){        // mass of primary was 0 or p = primary.  Return same particle without doing anything.
        return *p;
    }
    const double dt = sim->dt_last_done;
    double tau_a = INFINITY;
    double tau_e = INFINITY;
    double tau_inc = INFINITY;
    double tau_omega = INFINITY;
    double tau_Omega = INFINITY;
    
    rebx_get_param_check("tau_a", &tau_a, REBX_TYPE_DOUBLE, 1, p);
    rebx_get_param_check("tau_e", &tau_e, REBX_TYPE_DOUBLE, 1, p);
    rebx_get_param_check("tau_inc", &tau_inc, REBX_TYPE_DOUBLE, 1, p);
    rebx_get_param_check("tau_omega", &tau_omega, REBX_TYPE_DOUBLE, 1, p);
    rebx_get_param("tau_Omega", &tau_Omega, REBX_TYPE_DOUBLE, 1, p);
    
    const double a0 = o.a;
    const double e0 = o.e;
    const double inc0 = o.inc;

    o.a += a0*dt/tau_a;
    o.e += e0*dt/tau_e;
    o.inc += inc0*dt/tau_inc;
    o.omega += 2.*M_PI*dt/tau_omega;
    o.Omega += 2.*M_PI*dt/tau_Omega;
   
    if(tau_e < INFINITY){
        double p_param = 0;
        rebx_get_param("p", &p_param, REBX_TYPE_DOUBLE, 1, effect);
        o.a += 2.*a0*e0*e0*p_param*dt/tau_e; // Coupling term between e and a
    }
    return reb_tools_orbit_to_particle(sim->G, *primary, p->m, o.a, o.e, o.inc, o.Omega, o.omega, o.f);
}

void rebx_modify_orbits_direct(struct reb_simulation* const sim, struct rebx_effect* const effect){
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI;
    rebx_get_param("coordinates", &coordinates, REBX_TYPE_INT, 1, effect);
    
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebxtools_com_ptm(sim, effect, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_modify_orbits_direct);
}
