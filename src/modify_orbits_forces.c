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
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"
#include "core.h"

static struct reb_vec3d rebx_calculate_modify_orbits_forces(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* p, struct reb_particle* source){
    double tau_a = INFINITY;
    double tau_e = INFINITY;
    double tau_inc = INFINITY;
    
    const double* const tau_a_ptr = rebx_get_param_check(p, "tau_a", REBX_TYPE_DOUBLE);
    const double* const tau_e_ptr = rebx_get_param_check(p, "tau_e", REBX_TYPE_DOUBLE);
    const double* const tau_inc_ptr = rebx_get_param_check(p, "tau_inc", REBX_TYPE_DOUBLE);

    const double dvx = p->vx - source->vx;
    const double dvy = p->vy - source->vy;
    const double dvz = p->vz - source->vz;
    const double dx = p->x-source->x;
    const double dy = p->y-source->y;
    const double dz = p->z-source->z;
    const double r2 = dx*dx + dy*dy + dz*dz;
    
    if(tau_a_ptr != NULL){
        tau_a = *tau_a_ptr;
    }
    if(tau_e_ptr != NULL){
        tau_e = *tau_e_ptr;
    }
    if(tau_inc_ptr != NULL){
        tau_inc = *tau_inc_ptr;
    }
    
    /*const int* const timescales_in_orbital_periods = rebx_get_param_check(effect, "timescales_in_orbital_periods", REBX_TYPE_INT);
    
    if (timescales_in_orbital_periods != NULL){
        const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
        const double mu = sim->G*source->m;
        const double vcirc2 = mu/sqrt(r2);
        const double a = mu/(2.*vcirc2 - v2);
        if(a > 0){
            const double P = sqrt(a*a*a/mu);
            tau_a *= P;
            tau_e *= P;
            tau_inc *= P;
            const double* const dt_in_orbital_periods = rebx_get_param_check(p, "dt_in_orbital_periods", REBX_TYPE_DOUBLE);
            if(dt_in_orbital_periods != NULL){
                sim->dt = *dt_in_orbital_periods*P;
            }
        }
    }*/
    
    struct reb_vec3d a = {0};

    a.x =  dvx/(2.*tau_a);
    a.y =  dvy/(2.*tau_a);
    a.z =  dvz/(2.*tau_a);

    if (tau_e < INFINITY || tau_inc < INFINITY){
        const double vdotr = dx*dvx + dy*dvy + dz*dvz;
        const double prefac = 2*vdotr/r2/tau_e;
        a.x += prefac*dx;
        a.y += prefac*dy;
        a.z += prefac*dz + 2.*dvz/tau_inc;
    }
    return a;
}

void rebx_modify_orbits_forces(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* const particles, const int N){
    int* ptr = rebx_get_param_check(effect, "coordinates", REBX_TYPE_INT);
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI; // Default
    if (ptr != NULL){
        coordinates = *ptr;
    }
    
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    
    double* Edissipated_WHFAST = rebx_get_param_check(effect, "Edissipated_WHFAST", REBX_TYPE_DOUBLE);
    if(Edissipated_WHFAST != NULL){
        const double E0 = reb_tools_energy(sim);
        struct reb_particle* const ps = malloc(N*sizeof(*ps));
        memcpy(ps, particles, N*sizeof(*ps));
        rebx_reset_accelerations(sim->particles, N);
        rebx_com_force(sim, effect, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_modify_orbits_forces, particles, N);
        for(int i=0; i<N; i++){
            particles[i].vx += sim->dt*particles[i].ax;
            particles[i].vy += sim->dt*particles[i].ay;
            particles[i].vz += sim->dt*particles[i].az;
        }
        const double Ef = reb_tools_energy(sim);
        memcpy(particles, ps, N*sizeof(*ps));
        *Edissipated_WHFAST += Ef-E0;
        free(ps);
    }
    
    rebx_com_force(sim, effect, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_modify_orbits_forces, particles, N);
}
