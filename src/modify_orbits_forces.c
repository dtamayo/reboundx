/**
 * @file    modify_orbits_forces.c
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
 * Implementation Paper    *In progress*
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
 * =========================== ==================================================================
 * Field (C type)              Description
 * =========================== ==================================================================
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
#include "modify_orbits_forces.h"
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

static struct reb_vec3d rebx_calculate_modify_orbits_forces(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* p, struct reb_particle* source){
    double tau_a = INFINITY;
    double tau_e = INFINITY;
    double tau_inc = INFINITY;

    double* tau_a_ptr = rebx_get_param_double(p, "tau_a");
    double* tau_e_ptr = rebx_get_param_double(p, "tau_e");
    double* tau_inc_ptr = rebx_get_param_double(p, "tau_inc");

    if(tau_a_ptr != NULL){
        tau_a = *tau_a_ptr;
    }
    if(tau_e_ptr != NULL){
        tau_e = *tau_e_ptr;
    }
    if(tau_inc_ptr != NULL){
        tau_inc = *tau_inc_ptr;
    }

    const double dvx = p->vx - source->vx;
    const double dvy = p->vy - source->vy;
    const double dvz = p->vz - source->vz;

    struct reb_vec3d a = {0};

    a.x =  dvx/(2.*tau_a);
    a.y =  dvy/(2.*tau_a);
    a.z =  dvz/(2.*tau_a);

    if (tau_e < INFINITY || tau_inc < INFINITY){   // need h and e vectors for both types
        const double dx = p->x-source->x;
        const double dy = p->y-source->y;
        const double dz = p->z-source->z;
        const double r = sqrt ( dx*dx + dy*dy + dz*dz );
        const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
        const double prefac = 2*vr/r;
        a.x += prefac*dx/tau_e;
        a.y += prefac*dy/tau_e;
        a.z += prefac*dz/tau_e + 2.*dvz/tau_inc;
    }
    return a;
}

void rebx_modify_orbits_forces(struct reb_simulation* const sim, struct rebx_effect* const effect){
    int* ptr = rebx_get_param_int(effect, "coordinates");
    enum REBX_COORDINATES coordinates;
    if (ptr == NULL){
        coordinates = REBX_JACOBI;                  // Default
    }
    else{
        coordinates = *ptr;
    }
    
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_com_force(sim, effect, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_modify_orbits_forces);
}
/*
void rebx_modify_orbits_forces(struct reb_simulation* const sim, struct rebx_effect* const effect){
    const int N_real = sim->N - sim->N_var;

    struct reb_particle com;
    
    const enum REBX_COORDINATES coordinates;
    const enum REBX_COORDINATES* const ptr = rebx_get_param_int(effect, "coordinates");

    if (ptr == NULL){
        coordinates = JACOBI;                       // Default
    }
    else{
        coordinates = *ptr;
    }

    switch(coordinates){
    case JACOBI:
        com = reb_get_com(sim);                     // We start with outermost particle, so start with COM and peel off particles
        break;
    case BARYCENTRIC:
        com = reb_get_com(sim);                     // COM of whole system
        break;
    case HELIOCENTRIC:
        com = sim->particles[0];                    // Use the central body as the reference
        break;
    default:
        fprintf(stderr, "coordinates in parameters for modify_orbits_forces are not supported.\n");
        exit(1);
    }

    struct reb_particle* p0 = &sim->particles[0];
    for(int i=N_real-1;i>0;--i){
        struct reb_particle* p = &sim->particles[i];
        if(coordinates == JACOBI){
            rebxtools_update_com_without_particle(&com, p);
        }
        double tau_a = rebx_get_param_double(p, "tau_a");
        double tau_e = rebx_get_param_double(p, "tau_e");
        double tau_inc = rebx_get_param_double(p, "tau_inc");

        if(isnan(tau_a)){
            rebx_set_param_double(p, "tau_a", INFINITY);
            tau_a = INFINITY;
        }
        if(isnan(tau_e)){
            rebx_set_param_double(p, "tau_e", INFINITY);
            tau_e = INFINITY;
        }
        if(isnan(tau_inc)){
            rebx_set_param_double(p, "tau_inc", INFINITY);
            tau_inc = INFINITY;
        }

        const double dvx = p->vx - com.vx;
        const double dvy = p->vy - com.vy;
        const double dvz = p->vz - com.vz;

        double ax = 0.;
        double ay = 0.;
        double az = 0.;

        ax =  dvx/(2.*tau_a);
        ay =  dvy/(2.*tau_a);
        az =  dvz/(2.*tau_a);

        if (tau_e < INFINITY || tau_inc < INFINITY){   // need h and e vectors for both types
            const double dx = p->x-com.x;
            const double dy = p->y-com.y;
            const double dz = p->z-com.z;
            const double r = sqrt ( dx*dx + dy*dy + dz*dz );
            const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
            const double prefac = 2*vr/r;
            ax += prefac*dx/tau_e;
            ay += prefac*dy/tau_e;
            az += prefac*dz/tau_e + 2.*dvz/tau_inc;
        }
        p->ax += ax;
        p->ay += ay;
        p->az += az;
    }
    reb_move_to_com(sim);
}
*/
