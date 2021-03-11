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
 * dedge (double)               Yes         The position of the inner disc edge in code units 
 * hedge (double)               Yes         The aspect ratio at the inner disc edge; the disc edge width
 * sd0 (double)                 Yes         Initial disc surface density, must be in code units
 * zscale (double)              Yes         The aspect ratio at 1 code unit from the star
 * alpha (double)               Yes         Surface density profile of the disc
 * beta (double)                Yes         The flaring index; 1 means disc is irradiated by only the stellar flux
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
 * tau_a_red (double)           No          The planet trap to stop the migration at the inner disc edge
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"


const double rebx_calculate_planet_trap_type_I_mig(const double r, const double dedge, const double hedge){
    double tau_a_red;

    if (r > dedge*(1. + hedge)){
        tau_a_red = 1.;
    }

    else if (dedge*(1. - hedge) < r){
        tau_a_red =  5.5 * cos( ((dedge * (1. + hedge) - r ) * 2. * M_PI) / (4. * hedge * dedge) ) - 4.5;
    }

    else {
        tau_a_red = -10.;
    }

    return tau_a_red;
}

/* Calculating the t_wave: damping timescale or orbital evolution timescale from Tanaka & Ward 2004 (two papers with this equation slightly 
differently expressed). h = aspect ratio, sma = semi-major axis, sd0 = initial disc surface density, 
alpha = alpha a constant, sd = surface denisty to be calculated at every r, ms = stellar mass, mp = planet mass*/

const double rebx_calculating_damping_timescale(const double G, const double sd0, const double r, const double alpha, const double ms, const double mp, const double sma, const double h2){
    double sd;
    double t_wave;
    
    sd = sd0*pow(r, -alpha);
    t_wave = (sqrt(ms*ms*ms)*h2*h2)/(mp*sd*sqrt(sma*G));

    return t_wave;
}

/* Calculating the eccentricity damping timescale, eh=e/h*/

const double rebx_calculating_eccentricity_damping_timescale(const double wave, const double eh, const double ih){
    double t_e;
    t_e =  (wave/0.780) * (1. - (0.14*eh*eh) + (0.06*eh*eh*eh) + (0.18*eh*ih*ih)) ;
    
    return t_e;
}

/* Calculating the damping timescale of the semi-major axis; dampened as the planet moves inward. ih = inc/h */

const double rebx_calculating_semi_major_axis_damping_timescale(const double wave, const double eh, const double ih, const double h2, const double alpha){
    double Pe;
    double t_a;

    Pe = (1. + pow(eh/2.25, 1.2) +  pow(eh/2.84, 6.)) / (1. - pow(eh/2.02, 4.));
    t_a = ((2.*wave)/(2.7 + 1.1*alpha)) * (1/h2) * (Pe + (Pe/fabs(Pe)) * ((0.070*ih) + (0.085*ih*ih*ih*ih) - (0.080*eh*ih*ih)));

    return t_a;
}

/* Calculating the inclination damping timescale */
const double rebx_calculating_inclination_damping_timescale(const double wave, const double eh, const double ih){
    double t_i;
    t_i = (wave/0.544) * (1 - (0.30*ih*ih) + (0.24*ih*ih*ih) + (0.14*eh*eh*ih));

    return t_i;
}
static struct reb_vec3d rebx_calculate_modify_orbits_with_type_I_migration(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p, struct reb_particle* source){
    double invtau_a;
    double tau_e;
    double tau_inc;

    /* Default values for the parameters in case the user forgets to define them when using this code */
    double beta = 0.0;
    double zscale = 0.01;
    double sd0 = 0.0;
    double alpha = 0.0;
    double dedge = 0.0;
    double hedge = 0.0;

    /* Parameters that should be changed/set in Python notebook or C outside of this */
    const double* const dedge_ptr = rebx_get_param(sim->extras, force->ap, "inner_disc_edge");
    const double* const hedge_ptr = rebx_get_param(sim->extras, force->ap, "disc_edge_width");
    const double* const beta_ptr = rebx_get_param(sim->extras, force->ap, "flaring_index");
    const double* const alpha_ptr = rebx_get_param(sim->extras, force->ap, "alpha");
    const double* const sd0_ptr = rebx_get_param(sim->extras, force->ap, "initial_surface_density");
    const double* const zscale_ptr = rebx_get_param(sim->extras, force->ap, "aspect_ratio");

    /* This is accessing the calculated semi-major axis, eccentricity and inclination via modify_orbits_direct where they are calculated and returned*/
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

/* Setting some if statements for all the free parameters, so that if the user does not set a value the parameter will use the
default value given above, which is not correct (all zero) so the code will not work and the user will be able to realize why. These
parameters must be set in order for the code to work properly */
    if (dedge_ptr != NULL){
        dedge = *dedge_ptr;
    }
    if (hedge_ptr != NULL){
        hedge = *hedge_ptr;
    }
    if (beta_ptr != NULL){
        beta = *beta_ptr;
    }
    if (alpha_ptr != NULL){
        alpha = *alpha_ptr;
    }
    if (sd0_ptr != NULL){
        sd0 = *sd0_ptr;
    }
    if (zscale_ptr != NULL){
        zscale = *zscale_ptr;
    }

/* Calculating the aspect ratio evaluated at the position of the planet, r, with a default value of 0.01 which is
given/found when beta=0 at r = 1 code unit which is 1AU */

    const double h = (zscale) * pow(r2, beta/2); //In code units where r is in AU and beta can be the standard 0.25
    const double h2 = h*h;

/* Defining other variables needed in the different equations for the damping timescales */
    const double eh = e0/h;
    const double ih = inc0/h;

    const double G = sim->G;
    const double wave = rebx_calculating_damping_timescale(G, sd0, sqrt(r2), alpha, ms, mp, a0, h2);

    invtau_a = rebx_calculate_planet_trap_type_I_mig(a0, dedge, hedge)/(rebx_calculating_semi_major_axis_damping_timescale(wave, eh, ih, h2, alpha));
    tau_e = rebx_calculating_eccentricity_damping_timescale(wave, eh, ih);
    tau_inc = rebx_calculating_inclination_damping_timescale(wave, eh, ih);

    struct reb_vec3d a = {0};

    if (invtau_a != 0.0){
        a.x = -dvx*(invtau_a);
        a.y = -dvy*(invtau_a);
        a.z = -dvz*(invtau_a);
    }

    if (tau_e < INFINITY || tau_inc < INFINITY){
        const double vdotr = dx*dvx + dy*dvy + dz*dvz;
        const double prefac = -2*vdotr/r2/tau_e;
        a.x += prefac*dx;
        a.y += prefac*dy;
        a.z += prefac*dz - 2.*dvz/tau_inc;
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
