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
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"


const double rebx_calculate_planet_trap(const double r, const double h, const double dedge){
    double tau_a_red;

    if (r > dedge*(1.0 + h)){
        tau_a_red = 1.0;
    }

    else if (dedge*(1.0 - h) < r){
        tau_a_red =  5.5 * cos( ((dedge * (1.0 + h) - r ) * 2 * M_PI) / (4 * h * dedge) ) - 4.5;
    }

    else {
        tau_a_red = -10.0;
    }

    return tau_a_red;
}

/* Calculating the aspect ratio as a varying parameter, evaluated at the position of the planet, r, with a default value of 0.02 which is
given/found when beta=0. Ask Antoine how to work with default on C, should I do an if statement; if beta=0 or not given take Hr =0.02?? */

const double rebx_calculating_the_aspect_ratio(const double r, cosnt double beta){
    double Hr;
    Hr = 0.02*(3**(-beta))*(r**beta);
}

/* Calculating the t_wave: damping timescale or orbital evolution timescale from Tanaka & Ward 2004 (two papers with this equation slightly 
differently expressed) */
const double rebx_calculating_damping_timescale(...){
    double t_wave;
    t_wave = ((m_s*m_s)/(m*rebx_calculating_surface_density(...)*a0**2)) * (Hr**4) * rebx_calculating_the_angular_velocity(...);
}

/* Calculating the eccentricity damping timescale, all based on t_wave, the overall migration damping timescale*/
const double rebx_calculating_eccentricity_damping_timescale(...){
    double t_e;
    t_e =  (rebx_calculating_damping_timescale(...)/0.780) * (1.0 - 0.14*((e/Hr)**2) + 0.06*((e/Hr)**3));
}

/* Calculating the P(e) factor that will reverse the torque for high */
const double rebx_calculating_Pe_factor(const double Hr){
    double Pe;
    Pe = (1 + (e0/(2.25*Hr))**(1.2) + (e0/(2.84*Hr))**6) / (1 - (e0/(2.02*Hr))**4);
} 

/* Calculating the damping timescale of the semi-major axis, how it is dampened as the planet moves inward */
const double reb_calculating_semi_major_axis_damping_timescale(const double alpha, const double Hr){
    double t_a;
    t_a = ((2*rebx_calculating_damping_timescale(...))/(2.7 + 1.1*alpha)) * (Hr**2) * rebx_calculating_Pe_factor(...);
}
  

static struct reb_vec3d rebx_calculate_modify_orbits_with_type_I_migration(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p, struct reb_particle* source){
    double invtau_a = 0.0;
    double tau_e = INFINITY;
    double tau_inc = INFINITY;

    /* Parameters that can be changed/set in Python notebook or C outside of this */
    const double* const tau_a_ptr = rebx_get_param(sim->extras, p->ap, "tau_a");
    const double* const tau_e_ptr = rebx_get_param(sim->extras, p->ap, "tau_e");
    const double* const tau_inc_ptr = rebx_get_param(sim->extras, p->ap, "tau_inc");
    const double* const dedge = rebx_get_param(sim->extras, force->ap, "inner_disc_edge");
    const double* const h = rebx_get_param(sim->extras, force->ap, "disc_edge_width");
    const double* const Hr = rebx_get_param(sim->extras, force->ap, "aspect ratio");
    const double* const beta = rebx_get_param(sim->extras, force->ap, "beta");
    const double* const alpha = rebx_get_param(sim->extras, force->ap, "alpha");

    /* This is accessing the calculated semi-major axis, eccentricity and inclination via modify_orbits_direct where they are calculated and returned*/
    struct rebx_extras* const rebx = sim->extras;
    int err=0;
    struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *p, *primary, &err);
        
    const double a0 = o.a;
    const double e0 = o.e;
    const double inc0 = o.inc;


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
