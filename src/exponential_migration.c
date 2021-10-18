/**
 * @file    exponential_migration.c
 * @brief   Continuous velocity kicks leading to exponential change in the object's semimajor axis.
 * @author  Mohamad Ali-Dib <mma9132@nyu.edu>
 * 
 * @section     LICENSE
 * Copyright (c) 2021 Mohamad Ali-Dib
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
 * Author                   Mohamad Ali-Dib
 * Implementation Paper    `Ali-Dib et al., 2021 AJ <https://arxiv.org/abs/2104.04271>`_.
 * Based on                `Hahn & Malhotra 2005 <https://ui.adsabs.harvard.edu/abs/2005AJ....130.2392H/abstract>`_.
 * C Example               :ref:`c_example_exponential_migration`
 * Python Example          `ExponentialMigration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/ExponentialMigration.ipynb>`_.
 * ======================= ===============================================
 * 
 * Continuous velocity kicks leading to exponential change in the object's semimajor axis. 
 * One of the standard prescriptions often used in Neptune migration & Kuiper Belt formation models.
 * Does not directly affect the eccentricity or inclination of the object.
 * 
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * em_tau_a (double)              Yes          Semimajor axis exponential growth/damping timescale
 * em_aini (double)               Yes          Object's initial semimajor axis
 * em_afin (double)               Yes          Object's final semimajor axis
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

static struct  reb_vec3d rebx_calculate_modify_orbits_forces_new(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p,  struct reb_particle* source){

   struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, *p, *source);

    double em_tau_a = INFINITY;
    double em_aini = 24.;
    double em_afin = 30.;    

    const double* const em_tau_a_ptr = rebx_get_param(sim->extras, p->ap, "em_tau_a");
    const double* const em_ainipoint = rebx_get_param(sim->extras, p->ap, "em_aini");
    const double* const em_afinpoint = rebx_get_param(sim->extras, p->ap, "em_afin");

    const double dvx = p->vx - source->vx;
    const double dvy = p->vy - source->vy;
    const double dvz = p->vz - source->vz;
    const double dx = p->x-source->x;
    const double dy = p->y-source->y;
    const double dz = p->z-source->z;
    const double r2 = dx*dx + dy*dy + dz*dz;
    
    if(em_tau_a_ptr != NULL){
        em_tau_a = *em_tau_a_ptr;
    }
    if(em_ainipoint != NULL){
        em_aini = *em_ainipoint;
    }
    if(em_afinpoint != NULL){
        em_afin = *em_afinpoint;
    }
    
    
    struct reb_vec3d a = {0};

    a.x =  (dvx/(2.*em_tau_a))*((em_afin - em_aini)/(o.a))*exp(-(sim->t) / em_tau_a);
    a.y =  (dvy/(2.*em_tau_a))*((em_afin - em_aini)/(o.a))*exp(-(sim->t) / em_tau_a);
    a.z =  (dvz/(2.*em_tau_a))*((em_afin - em_aini)/(o.a))*exp(-(sim->t) / em_tau_a);


    return a;
}


void rebx_exponential_migration(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    int* ptr = rebx_get_param(sim->extras, force->ap, "coordinates");
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI; // Default
    if (ptr != NULL){
        coordinates = *ptr;
    }
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_com_force(sim, force, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_modify_orbits_forces_new, particles, N);
}
