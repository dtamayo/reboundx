/**
 * @file    type_I_migration.c
 * @brief   Type I migration 
 * @author  Kaltrina Kajtazi <1kaltrinakajtazi@gmail.com>, Gabriele Pichierri <gabrielepichierri@gmail.com>
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
 * $Orbit Modifications$       // Effect category 
 * 
 * ======================= ===============================================
 * Authors                 Kajtazi, Kaltrina and D. Petit, C. Antoine
 * Implementation Paper    `Kajtazi et al 2022 <https://ui.adsabs.harvard.edu/abs/2022arXiv221106181K/abstract>`_.
 * Based on                `Cresswell & Nelson 2008 <https://ui.adsabs.harvard.edu/abs/2008A%26A...482..677C/abstract>`_, and `Pichierri et al 2018 <https://ui.adsabs.harvard.edu/abs/2018CeMDA.130...54P/abstract>`_.
 * C example               :ref:`c_example_type_I_migration`
 * Python example          `TypeIMigration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TypeIMigration.ipynb>`_.
 * ======================= ===============================================
 * 
 * This applies Type I migration, damping eccentricity, angular momentum and inclination.
 * The base of the code is the same as the modified orbital forces one written by D. Tamayo, H. Rein.
 * It also allows for parameters describing an inner disc edge, modeled using the implementation in inner_disk_edge.c.
 * Note that this code is not machine independent since power laws were not possible to avoid all together.
 *
 * **Effect Parameters**
 * 
 * ===================================== =========== ==================================================================================================================
 * Field (C type)                        Required    Description
 * ===================================== =========== ==================================================================================================================
 * ide_position (double)                 No          The position of the inner disk edge in code units 
 * ide_width (double)                    No          The disk edge width (planet will stop within ide_width of ide_position)
 * tIm_surface_density_1 (double)        Yes         Disk surface density at one code unit from the star; used to find the surface density at any distance from the star
 * tIm_scale_height_1 (double)           Yes         The scale height at one code unit from the star; used to find the aspect ratio at any distance from the star
 * tIm_surface_density_exponent (double) Yes         Exponent of disk surface density, indicative of the surface density profile of the disk
 * tIm_flaring_index (double)            Yes         The flaring index; 1 means disk is irradiated by only the stellar flux
 * ===================================== =========== ==================================================================================================================
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"


/* Calculating the t_wave: damping timescale or orbital evolution timescale, from Tanaka & Ward 2004. 
h = aspect ratio, h2 = aspect ratio squared, sma = semi-major axis, sd = disk surface denisty to be calculated at every r, ms = stellar mass, mp = planet mass */

const double rebx_calculate_damping_timescale(const double G, const double sd0, const double r, const double s, const double ms, const double mp, const double sma, const double h2){
    double sd;
    double t_wave;
    
    sd = sd0*pow(r, -s);
    t_wave = (sqrt(ms*ms*ms)*h2*h2)/(mp*sd*sqrt(sma*G));

    return t_wave;
}

/* Calculating the eccentricity damping timescale t_e = -e/(de/dt), from Cresswell & Nelson 2008. 
eh=e/h, ih = i/h, wave is a funcion variable name which will be the t_wave function*/

const double rebx_calculate_eccentricity_damping_timescale(const double wave, const double eh, const double ih){
    double t_e;
    t_e =  (wave/0.780) * (1. - (0.14*eh*eh) + (0.06*eh*eh*eh) + (0.18*eh*ih*ih)) ;
    
    return t_e;
}

/* Calculating the migration timescale t_mig = - angmom/torque, from Cresswell & Nelson 2008*/

const double rebx_calculate_migration_timescale(const double wave, const double eh, const double ih, const double h2, const double s){
    double Pe;
    double t_mig;
    double term;
    double term2;
    double term3;
    term = (eh/2.25);
    term2 = (eh/2.84) * (eh/2.84);
    term3 = (eh/2.02) * (eh/2.02);
    Pe = (1. + pow(term, 1.2) +  term2*term2*term2) / (1. - term3*term3);
    t_mig = ((2.*wave)/(2.7 + 1.1*s)) * (1/h2) * (Pe + (Pe/fabs(Pe)) * ((0.070*ih) + (0.085*ih*ih*ih*ih) - (0.080*eh*ih*ih)));

    return t_mig;
}

/* Calculating the inclination damping timescale t_i = -i/(di/dt), from Cresswell & Nelson 2008*/

const double rebx_calculate_inclination_damping_timescale(const double wave, const double eh, const double ih){
    double t_i;
    t_i = (wave/0.544) * (1 - (0.30*ih*ih) + (0.24*ih*ih*ih) + (0.14*eh*eh*ih));

    return t_i;
}

static struct reb_vec3d rebx_calculate_modify_orbits_with_type_I_migration(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p, struct reb_particle* source){
    double invtau_mig;
    double tau_e;
    double tau_inc;

    /* Default values for the parameters in case the user forgets to define them when using this code */
    double beta = 0.0;
    double h0 = 0.01;
    double sd0 = 0.0;
    double s = 0.0;
    double dedge = 0.0;
    double hedge = 0.0;

    /* Parameters that should be changed/set in Python notebook or in C outside of this */
    const double* const dedge_ptr = rebx_get_param(sim->extras, force->ap, "ide_position");
    const double* const hedge_ptr = rebx_get_param(sim->extras, force->ap, "ide_width");
    const double* const beta_ptr = rebx_get_param(sim->extras, force->ap, "tIm_flaring_index");
    const double* const s_ptr = rebx_get_param(sim->extras, force->ap, "tIm_surface_density_exponent");
    const double* const sd0_ptr = rebx_get_param(sim->extras, force->ap, "tIm_surface_density_1");
    const double* const h0_ptr = rebx_get_param(sim->extras, force->ap, "tIm_scale_height_1");

    /* Accessing the calculated semi-major axis, eccentricity and inclination for each integration step, via modify_orbits_direct where they are calculated and returned*/
    int err=0;
    struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *p, *source, &err);
  
    const double a0 = o.a;
    const double e0 = o.e;
    const double inc0 = o.inc;
    const double mp = p->m;  
    const double ms = source->m;

    const double dvx = p->vx - source->vx;
    const double dvy = p->vy - source->vy;
    const double dvz = p->vz - source->vz;
    const double dx = p->x-source->x;
    const double dy = p->y-source->y;
    const double dz = p->z-source->z;
    const double r2 = dx*dx + dy*dy + dz*dz;

    if (beta_ptr != NULL){
        beta = *beta_ptr;
    }
    if (s_ptr != NULL){
        s = *s_ptr;
    }
    if (sd0_ptr != NULL){
        sd0 = *sd0_ptr;
    }
    if (h0_ptr != NULL){
        h0 = *h0_ptr;
    }
    if (dedge_ptr != NULL){
        dedge = *dedge_ptr;
    }
    if (hedge_ptr != NULL){
        hedge = *hedge_ptr;
    }

    /* Calculating the aspect ratio evaluated at the position of the planet, r and defining other variables */

    const double h = (h0) * pow(r2, beta/2); 
    const double h2 = h*h;

    const double eh = e0/h;
    const double ih = inc0/h;

    const double G = sim->G;
    const double wave = rebx_calculate_damping_timescale(G, sd0, sqrt(r2), s, ms, mp, a0, h2);
    invtau_mig = rebx_calculate_planet_trap(a0, dedge, hedge)/(rebx_calculate_migration_timescale(wave, eh, ih, h2, s));
    tau_e = rebx_calculate_eccentricity_damping_timescale(wave, eh, ih);
    tau_inc = rebx_calculate_inclination_damping_timescale(wave, eh, ih);

    struct reb_vec3d a = {0};

    if (invtau_mig != 0.0){
        a.x = -dvx*(invtau_mig);
        a.y = -dvy*(invtau_mig);
        a.z = -dvz*(invtau_mig);
    }

    if (tau_e < INFINITY || tau_inc < INFINITY){
        const double vdotr = dx*dvx + dy*dvy + dz*dvz;
        const double prefac = -2*vdotr/r2/tau_e;
        a.x += prefac*dx;
        a.y += prefac*dy;
        a.z += prefac*dz - 2*dvz/tau_inc;
    }
    return a;
}

void rebx_modify_orbits_with_type_I_migration(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    int* ptr = rebx_get_param(sim->extras, force->ap, "coordinates");
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI; // Default
    if (ptr != NULL){
        coordinates = *ptr;
    }
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_com_force(sim, force, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_modify_orbits_with_type_I_migration, particles, N);
}
