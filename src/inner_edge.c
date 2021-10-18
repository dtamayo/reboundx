/**
 * @file    inner_disk_edge.c
 * @brief   Apply an inner disk edge at a chosen location, when integating planetary motion
 * @author  Kaltrina Kajtazi <1kaltrinakajtazi@gmail.com>
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
 
 * $Inner disc edge$       // Effect category 
 * 
 * ======================= ============================================================================================
 * Authors                 Kajtazi, Kaltrina and D. Petit, C. Antoine
 * Implementation Paper    `Kajtazi et al. in prep.
 * Based on                `Pichierri et al 2018 <https://ui.adsabs.harvard.edu/abs/2018CeMDA.130...54P/abstract>.
 * ======================= ============================================================================================
 * 
 * This applies an inner disc edge that functions as a planet trap. Within its width the planet's migration is reversed 
 * by an opposite and roughly equal magnitude torque. Thus, stopping further migration and trapping the planet within 
 * the width of the trap. The base used here is modified_orbital_forces script written by D. Tamayo, H. Rein.
 * 
 * This implementation should work with any migration/effect not just Type I migration or constant migration. 
 * Other precriptions have not been tested but should work fine, as long as that migration prescription can be given 
 * in terms of the timescales of change in orbital elements and applied through accelerations as done here. 
 * However, the code is not machine idependent due to the use of power laws, which cannot be avoided altogether in this case. 
 * 
 * **Effect Parameters**
 * 
 * ============================ =========== ===================================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ===================================================================================
 * dedge (double)               Yes         The position of the inner disk edge in code units 
 * hedge (double)               Yes         The aspect ratio at the inner disk edge; the disk edge width
 * ============================ =========== ===================================================================================
 *
 * **Particle Parameters**
 *
 * One can pick and choose which particles have which parameters set.  
 *
 * ============================ =========== ===================================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ===================================================================================
 * tau_a (double)               No          Semimajor axis exponential growth/damping timescale
 * tau_e (double)               No          Eccentricity exponential growth/damping timescale
 * tau_inc (double)             No          Inclination axis exponential growth/damping timescale
 * tau_a_red (double)           No          Planet trap function to stop further migration once the inner disc edge is reached
 * ============================ =========== ===================================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

/* A planet trap that is active only at the inner disc edge, to reverse the planetary migration and prevent migration onto the star */
const double rebx_calculating_planet_trap(const double r, const double hedge, const double dedge){
    double tau_a_red;

    if (r > dedge*(1.0 + hedge)){
        tau_a_red = 1.0;
    }

    else if (dedge*(1.0 - hedge) < r){
        tau_a_red =  5.5 * cos( ((dedge * (1.0 + hedge) - r ) * 2 * M_PI) / (4 * hedge * dedge) ) - 4.5;
    }

    else {
        tau_a_red = -10.0;
    }

    return tau_a_red;
}

static struct reb_vec3d rebx_calculate_modify_orbits_with_inner_disc_edge(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p, struct reb_particle* source){
    double invtau_a = 0.0;
    double tau_e = INFINITY;
    double tau_inc = INFINITY;
    
    const double* const tau_a_ptr = rebx_get_param(sim->extras, p->ap, "tau_a");
    const double* const tau_e_ptr = rebx_get_param(sim->extras, p->ap, "tau_e");
    const double* const tau_inc_ptr = rebx_get_param(sim->extras, p->ap, "tau_inc");
    const double* const dedge = rebx_get_param(sim->extras, force->ap, "inner_disc_edge");
    const double* const hedge = rebx_get_param(sim->extras, force->ap, "disc_edge_width");

    const double dvx = p->vx - source->vx;
    const double dvy = p->vy - source->vy;
    const double dvz = p->vz - source->vz;
    const double dx = p->x-source->x;
    const double dy = p->y-source->y;
    const double dz = p->z-source->z;
    const double r2 = dx*dx + dy*dy + dz*dz;

    int err=0;
    struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *p, *source, &err);

    const double a0 = o.a;

    if(tau_a_ptr != NULL){
        invtau_a = rebx_calculating_planet_trap(a0, *hedge, *dedge)/(*tau_a_ptr);
    }
    if(tau_e_ptr != NULL){
        tau_e = *tau_e_ptr;
    }
    if(tau_inc_ptr != NULL){
        tau_inc = *tau_inc_ptr;
    }
    
    struct reb_vec3d a = {0};

    a.x =  dvx*invtau_a/(2.);
    a.y =  dvy*invtau_a/(2.);
    a.z =  dvz*invtau_a/(2.);

    if (tau_e < INFINITY || tau_inc < INFINITY){
        const double vdotr = dx*dvx + dy*dvy + dz*dvz;
        const double prefac = 2*vdotr/r2/tau_e;
        a.x += prefac*dx;
        a.y += prefac*dy;
        a.z += prefac*dz + 2.*dvz/tau_inc;
    }
    return a;
}


void rebx_modify_orbits_with_inner_disc_edge(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    int* ptr = rebx_get_param(sim->extras, force->ap, "coordinates");
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI; // Default
    if (ptr != NULL){
        coordinates = *ptr;
    }
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_com_force(sim, force, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_modify_orbits_with_inner_disc_edge, particles, N);
}

