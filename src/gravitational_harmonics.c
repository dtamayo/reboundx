/** * @file gravitational_harmonics.c
 * @brief   Add J2n gravitational harmonics to particles
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
 * $Gravity Fields$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 M. Broz
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_. 
 * Based on                None
 * C Example               :ref:`c_example_J2`
 * Python Example          `J2.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/J2.ipynb>`_.
 * ======================= ===============================================
 * 
 * Adds azimuthally symmetric gravitational harmonics (J2, J4) to bodies in the simulation.
 *
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * J2 (double)                  No          J2 coefficient
 * J4 (double)                  No          J4 coefficient
 * R_eq (double)                No          Equatorial radius of nonspherical body used for calculating Jn harmonics
 * Omega (reb_vec3d)            No          Angular rotation frequency (Omega_x, Omega_y, Omega_z)
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_j2(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){

    const double G = sim->G;
    struct rebx_extras* const rebx = sim->extras;

    for (int i=0; i<N; i++){
        const double* const J2 = rebx_get_param(rebx, particles[i].ap, "J2");
        if (J2 == NULL){
            continue;
        }
        const double* const R_eq = rebx_get_param(rebx, particles[i].ap, "R_eq");
        if (R_eq == NULL){
            continue;
        }
        const struct reb_vec3d* Omega = rebx_get_param(rebx, particles[i].ap, "Omega");
        if (Omega == NULL){
            continue;
        }
        const struct reb_particle pi = particles[i];

        const double omega2 = Omega->x*Omega->x + Omega->y*Omega->y + Omega->z*Omega->z;
        const double omega = sqrt(omega2);
        struct reb_vec3d s = {0};
        s.x = Omega->x/omega;
        s.y = Omega->y/omega;
        s.z = Omega->z/omega;

        /* new coordinate basis (body-fixed) */
        struct reb_vec3d hatu = {0};
        struct reb_vec3d hatv = {0};
        struct reb_vec3d hatw = {0};
        hatw.x = s.x;
        hatw.y = s.y;
        hatw.z = s.z;
        double fac = sqrt(s.x*s.x + s.y*s.y);
        if (fac != 0.0) {
            hatu.x = -s.y/fac;
            hatu.y = s.x/fac;
            hatu.z = 0.0;
        } else {
            hatu.x = 1.0;
            hatu.y = 0.0;
            hatu.z = 0.0;
        }
        hatv.x = -(hatu.y*hatw.z - hatu.z*hatw.y);
        hatv.y = -(hatu.z*hatw.x - hatu.x*hatw.z);
        hatv.z = -(hatu.x*hatw.y - hatu.y*hatw.x);

        /* old basis (') in new coordinates */
        struct reb_vec3d hatx_ = {0};
        struct reb_vec3d haty_ = {0};
        struct reb_vec3d hatz_ = {0};
        hatx_.x = hatu.x;
        hatx_.y = hatv.x;
        hatx_.z = hatw.x;
        haty_.x = hatu.y;
        haty_.y = hatv.y;
        haty_.z = hatw.y;
        hatz_.x = hatu.z;
        hatz_.y = hatv.z;
        hatz_.z = hatw.z;

        for (int j=0; j<N; j++){
            if (j == i){
                continue;
            }
            const struct reb_particle pj = particles[j];

            const double dx = pj.x - pi.x;
            const double dy = pj.y - pi.y;
            const double dz = pj.z - pi.z;
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double r = sqrt(r2);

            /* new coordinates */
            const double du = hatu.x*dx + hatu.y*dy + hatu.z*dz;
            const double dv = hatv.x*dx + hatv.y*dy + hatv.z*dz;
            const double dw = hatw.x*dx + hatw.y*dy + hatw.z*dz;
            const double costheta = dw/r;

            const double f1 = 3.0/2.0*G*pi.m*(*J2)*(*R_eq)*(*R_eq)/r2/r2/r;
            const double f2 = 5.0*costheta*costheta - 1.0;
            const double f3 = f2 - 2.0;

            const double au = f1*f2*du;
            const double av = f1*f2*dv;
            const double aw = f1*f3*dw;

            /* old coordinates */
            const double ax = hatx_.x*au + haty_.x*av + hatz_.x*aw;
            const double ay = hatx_.y*au + haty_.y*av + hatz_.y*aw;
            const double az = hatx_.z*au + haty_.z*av + hatz_.z*aw;

            particles[j].ax += ax;
            particles[j].ay += ay;
            particles[j].az += az;

            fac = pj.m/pi.m;

            particles[i].ax -= fac*ax;
            particles[i].ay -= fac*ay;
            particles[i].az -= fac*az;
        }
    }
}

void rebx_gravitational_harmonics(struct reb_simulation* const sim, struct rebx_force* const gh, struct reb_particle* const particles, const int N){
    rebx_j2(sim, gh, particles, N);
}


