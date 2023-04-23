/** * @file lense_thirring.c
 * @brief   Add Lense-Thirring effect to particles
 * @author  Arya Akmal <akmala@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2023 Arya Akmal
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
 * $General Relativity$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 A. Akmal
 * Implementation Paper    None
 * Based on                `Park et al. <https://iopscience.iop.org/article/10.3847/1538-3881/abd414/>`_.
 * C Example               :ref:`c_example_lense_thirring`
 * Python Example          `LT.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/LT.ipynb>`_.
 * ======================= ===============================================
 * 
 * Adds Lense-Thirring effect due to massive rotating bodies in the simulation.
 *
 * **Effect Parameters**
 * 
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * c (double)                   Yes         Speed of light in the units used for the simulation.
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * omega (double)               No          rotation rate
 * R_eq (double)                No          Equatorial radius of source body
 * C_fac (double)               No          Moment of Inertia of source body over MR^2
 * p_hat_x (double)             No          x-component of spin-pole unit vector
 * p_hat_y (double)             No          y-component of spin-pole unit vector
 * p_hat_z (double)             No          z-component of spin-pole unit vector
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_LT_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double omega, const double R_eq, const double C_fac, const double p_hat_x, const double p_hat_y, const double p_hat_z, const int source_index, const double C2){
    const struct reb_particle source = particles[source_index];
    const double G = sim->G;
    const double gamma = 1.000021;   //hard-coded Eddington-Robertson-Shiff parameter for now
    for (int i=0; i<N; i++){
        if(i == source_index){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
        const double r3 = r2*r;
        const double dvx = p.vx - source.vx;
        const double dvy = p.vy - source.vy;
        const double dvz = p.vz - source.vz;
        const double Jx = C_fac*source.m * R_eq*R_eq*omega*p_hat_x ;
        const double Jy = C_fac*source.m * R_eq*R_eq*omega*p_hat_y ;
        const double Jz = C_fac*source.m * R_eq*R_eq*omega*p_hat_z ;
        const double Omega_fac = (1.+gamma)*G/2/C2;
        const double Omega_x = Omega_fac*(-Jx +3.*(Jx*dx+Jy*dy+Jz*dz)*dx/r2)/r3;
        const double Omega_y = Omega_fac*(-Jy +3.*(Jx*dx+Jy*dy+Jz*dz)*dy/r2)/r3;
        const double Omega_z = Omega_fac*(-Jz +3.*(Jx*dx+Jy*dy+Jz*dz)*dz/r2)/r3;

        particles[i].ax += 2.*Omega_y*dvz - Omega_z*dvy;
        particles[i].ay += 2.*Omega_z*dvx - Omega_x*dvz;
        particles[i].az += 2.*Omega_x*dvy - Omega_y*dvx;
        particles[source_index].ax -= 2.*Omega_y*dvz - Omega_z*dvy; 
        particles[source_index].ay -= 2.*Omega_z*dvx - Omega_x*dvz;
        particles[source_index].az -= 2.*Omega_x*dvy - Omega_y*dvx;
    }
}

void rebx_lense_thirring(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    double* c = rebx_get_param(sim->extras, force->ap, "lt_c");
    if (c == NULL){
        reb_error(sim, "REBOUNDx Error: Need to set speed of light in LT effect.  See examples in documentation.\n");
    }
    const double C2 = (*c)*(*c);

    for (int i=0; i<N; i++){
        const double* omega = rebx_get_param(rebx, particles[i].ap, "lt_rot_rate");
        if(omega != NULL){
           const double* R_eq = rebx_get_param(rebx, particles[i].ap, "lt_R_eq");
           if (R_eq != NULL){
               const double* C_fac = rebx_get_param(rebx, particles[i].ap, "lt_Mom_I_fac");
               if(C_fac != NULL){
                  const double* p_hat_x = rebx_get_param(rebx, particles[i].ap, "lt_p_hatx");
                  if(p_hat_x != NULL){
                     const double* p_hat_y = rebx_get_param(rebx, particles[i].ap, "lt_p_haty");
                     if(p_hat_y != NULL){
                        const double* p_hat_z = rebx_get_param(rebx, particles[i].ap, "lt_p_hatz");
                        if(p_hat_z != NULL){
                           rebx_calculate_LT_force(sim, particles, N, *omega, *R_eq, *C_fac, *p_hat_x, *p_hat_y, *p_hat_z, i, C2);
                        }
                     }
                  }
               }
           }
        }
    }
}


