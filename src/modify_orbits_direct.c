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
 * =========================== ==================================================================
 * Field (C type)              Description
 * =========================== ==================================================================
 * p (double)                  Coupling parameter between eccentricity and semimajor axis evolution
 *                             (see Deck & Batygin 2015). `p=0` corresponds to no coupling, `p=1` to
 *                             eccentricity evolution at constant angular momentum.
 * coordinates (enum)          Type of elements to use for modification (Jacobi, barycentric or heliocentric).
 *                             See the examples for usage.
 * =========================== ==================================================================
 * 
 * **Particle Parameters**
 * 
 * One can pick and choose which particles have which parameters set.  
 * For each particle, any unset parameter is ignored.
 * 
 * =========================== ======================================================
 * Name (C type)               Description
 * =========================== ======================================================
 * tau_a (double)              Semimajor axis exponential growth/damping timescale
 * tau_e (double)              Eccentricity exponential growth/damping timescale
 * tau_inc (double)            Inclination axis exponential growth/damping timescale
 * tau_Omega (double)          Period of linear nodal precession/regression
 * tau_omega (double)          Period of linear apsidal precession/regression
 * =========================== ======================================================
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "modify_orbits_direct.h"
#include "rebound.h"
#include "reboundx.h"

#define TWOPI 6.2831853071795862

struct rebx_params_modify_orbits_direct* rebx_add_modify_orbits_direct(struct rebx_extras* rebx){
	struct rebx_params_modify_orbits_direct* params = malloc(sizeof(*params));
	params->p = 0;
    params->coordinates = JACOBI;
    rebx_add_post_timestep_modification(rebx, params, "modify_orbits_direct", rebx_modify_orbits_direct);
    return params;
}

void rebx_modify_orbits_direct(struct reb_simulation* const sim, struct rebx_effect* const effect){
    const struct rebx_params_modify_orbits_direct* const params = effect->paramsPtr;
    const int N_real = sim->N - sim->N_var;

    struct reb_particle com;
    switch(params->coordinates){
    case JACOBI:
        com = reb_get_com(sim);                             // We start with outermost particle, so start with COM and peel off each particle
        break;
    case BARYCENTRIC:
        com = reb_get_com(sim);                             // COM of whole system
        break;
    case HELIOCENTRIC:
        com = sim->particles[0];                            // Use the central body as the reference
        break;
    default:
        fprintf(stderr, "coordinates in parameters for modify_orbits_direct are not supported.\n");
        exit(1);
    }

    for(int i=N_real-1;i>0;--i){
        struct reb_particle* p = &(sim->particles[i]);
        if(params->coordinates == JACOBI){
            rebxtools_update_com_without_particle(&com, p);
        }
        struct reb_orbit o = rebxtools_particle_to_orbit(sim->G, p, &com);
        const double dt = sim->dt_last_done;
        const double tau_a = rebx_get_param_double(p,"tau_a");
        const double tau_e = rebx_get_param_double(p,"tau_e");
        const double tau_inc = rebx_get_param_double(p,"tau_inc");
        const double tau_omega = rebx_get_param_double(p,"tau_omega");
        const double tau_Omega = rebx_get_param_double(p,"tau_Omega");
        
        const double a0 = o.a;
        const double e0 = o.e;
        const double inc0 = o.inc;

        if(!isnan(tau_a)){
            o.a += a0*(dt/tau_a);
        }
        if(!isnan(tau_e)){
            o.e += e0*(dt/tau_e);
        }
        if(!isnan(tau_inc)){
            o.inc += inc0*(dt/tau_inc);
        }
        if(!isnan(tau_omega)){
            o.omega += TWOPI*dt/tau_omega;
        }
        if(!isnan(tau_Omega)){
            o.Omega += TWOPI*dt/tau_Omega;
        }
        if(!isnan(tau_e)){
            o.a += 2.*a0*e0*e0*params->p*dt/tau_e; // Coupling term between e and a
        }

        rebxtools_orbit2p(sim->G, p, &com, &o);
    }
}


