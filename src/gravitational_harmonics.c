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
 * C Example               :ref:`c_example_j2`
 * Python Example          `J2.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/J2.ipynb>`_.
 * ======================= ===============================================
 * 
 * Allows the user to add azimuthally symmetric gravitational harmonics (J2, J4) to bodies in the simulation. 
 * These interact with all other bodies in the simulation (treated as point masses). 
 * The implementation allows the user to specify an arbitrary spin axis orientation for each oblate body, which defines the axis of symmetry.
 * This is specified through the angular rotation rate vector Omega.
 * The rotation rate Omega is not currently used other than to specify the spin axis orientation.
 * In particular, the current implementation applies the appropriate torque from the body's oblateness to the orbits of all the other planets, but does not account for the equal and opposite torque on the body's spin angular momentum.
 * The bodies spins therefore remain constant in the current implementation.
 * This is a good approximation in the limit where the bodies' spin angular momenta are much greater than the orbital angular momenta involved.
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

#define DEFAULTOMEGA {0.0, 0.0, 1.0}

inline void j2_func(double G, double m, const double* J2, const double* R_eq, double r, double r2, double costheta2, double du, double dv, double dw, double* au, double* av, double* aw) {

    if (J2 == NULL) { return; }
    if (*J2 == 0.0) { return; }

    const double f1 = 3.0/2.0*G*m*(*J2)*(*R_eq)*(*R_eq)/r2/r2/r;
    const double f2 = 5.0*costheta2 - 1.0;
    const double f3 = f2 - 2.0;

    *au += f1*f2*du;
    *av += f1*f2*dv;
    *aw += f1*f3*dw;

    return;
}

inline void j4_func(double G, double m, const double* J4, const double* R_eq, double r, double r2, double costheta2, double du, double dv, double dw, double* au, double* av, double* aw) {

    if (J4 == NULL) { return; }
    if (*J4 == 0.0) { return; }

    const double f1 = 5.0/8.0*G*m*(*J4)*(*R_eq)*(*R_eq)*(*R_eq)*(*R_eq)/r2/r2/r2/r;
    const double f2 = 63.0*costheta2*costheta2 - 42.0*costheta2 + 3.0;
    const double f3 = f2 - 28.0*costheta2 + 12.0;

    *au += f1*f2*du;
    *av += f1*f2*dv;
    *aw += f1*f3*dw;

    return;
}

inline void uvw(struct reb_vec3d Omega, struct reb_vec3d* hatu, struct reb_vec3d* hatv, struct reb_vec3d* hatw) {

    const double omega2 = Omega.x*Omega.x + Omega.y*Omega.y + Omega.z*Omega.z;
    const double omega = sqrt(omega2);

    struct reb_vec3d s = {0};
    s.x = Omega.x/omega;
    s.y = Omega.y/omega;
    s.z = Omega.z/omega;

    hatw->x = s.x;
    hatw->y = s.y;
    hatw->z = s.z;

    double fac = sqrt(s.x*s.x + s.y*s.y);
    if (fac != 0.0) {
        hatu->x = -s.y/fac;
        hatu->y = s.x/fac;
        hatu->z = 0.0;
    } else {
        hatu->x = 1.0;
        hatu->y = 0.0;
        hatu->z = 0.0;
    }

    hatv->x = -(hatu->y*hatw->z - hatu->z*hatw->y);
    hatv->y = -(hatu->z*hatw->x - hatu->x*hatw->z);
    hatv->z = -(hatu->x*hatw->y - hatu->y*hatw->x);

    return;
}

void rebx_gravitational_harmonics(struct reb_simulation* const sim, struct rebx_force* const gh, struct reb_particle* const particles, const int N){
    const double G = sim->G;
    struct rebx_extras* const rebx = sim->extras;

    for (int i=0; i<N; i++){
        const double* const J2 = rebx_get_param(rebx, particles[i].ap, "J2");
        if (J2 == NULL){
            continue;
        }
        if (*J2 == 0.0){
            continue;
        }
        const double* const J4 = rebx_get_param(rebx, particles[i].ap, "J4");
        const double* const R_eq = rebx_get_param(rebx, particles[i].ap, "R_eq");
        if (R_eq == NULL){
            continue;
        }
        struct reb_vec3d Omega = DEFAULTOMEGA;
        const struct reb_vec3d* Omegaptr = rebx_get_param(rebx, particles[i].ap, "Omega");
        if (Omegaptr != NULL){
            Omega.x = Omegaptr->x;
            Omega.y = Omegaptr->y;
            Omega.z = Omegaptr->z;
        }
        const struct reb_particle pi = particles[i];

        /* new coordinate basis (body-fixed) */
        struct reb_vec3d hatu = {0};
        struct reb_vec3d hatv = {0};
        struct reb_vec3d hatw = {0};

        uvw(Omega, &hatu, &hatv, &hatw);

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
            const double costheta2 = costheta*costheta;

            double au = 0.0;
            double av = 0.0;
            double aw = 0.0;

            j2_func(G, pi.m, J2, R_eq, r, r2, costheta2, du, dv, dw, &au, &av, &aw);
            j4_func(G, pi.m, J4, R_eq, r, r2, costheta2, du, dv, dw, &au, &av, &aw);

            /* old coordinates */
            const double ax = hatx_.x*au + hatx_.y*av + hatx_.z*aw;
            const double ay = haty_.x*au + haty_.y*av + haty_.z*aw;
            const double az = hatz_.x*au + hatz_.y*av + hatz_.z*aw;

            particles[j].ax += ax;
            particles[j].ay += ay;
            particles[j].az += az;

            const double fac = pj.m/pi.m;

            particles[i].ax -= fac*ax;
            particles[i].ay -= fac*ay;
            particles[i].az -= fac*az;
        }
    }
}

inline void j2_potential_func(double G, double mi, double mj, const double* J2, const double* R_eq, double r, double r2, double costheta2, double* H) {

    if (J2 == NULL) { return; }
    if (*J2 == 0.0) { return; }

    const double f1 = G*mi*mj*(*J2)*(*R_eq)*(*R_eq)/r2/r;
    const double f2 = 1.0/2.0*(3.0*costheta2 - 1.0);

    *H += f1*f2;

    return;
}

inline void j4_potential_func(double G, double mi, double mj, const double* J4, const double* R_eq, double r, double r2, double costheta2, double* H) {

    if (J4 == NULL) { return; }
    if (*J4 == 0.0) { return; }

    const double f1 = G*mi*mj*(*J4)*(*R_eq)*(*R_eq)*(*R_eq)*(*R_eq)/r2/r2/r;
    const double f2 = 1.0/8.0*(35.0*costheta2*costheta2 - 30.0*costheta2 + 3.0);

    *H += f1*f2;

    return;
}

double rebx_gravitational_harmonics_potential(struct rebx_extras* const rebx){
    if (rebx->sim == NULL){
        rebx_error(rebx, "");
        return 0;
    }
    const struct reb_simulation* const sim = rebx->sim;
    const struct reb_particle* const particles = sim->particles;
    const double G = sim->G;
    const int N = sim->N - sim->N_var;
    double H = 0.0;

    for (int i=0; i<N; i++){
        const double* const J2 = rebx_get_param(rebx, particles[i].ap, "J2");
        if (J2 == NULL){
            continue;
        }
        if (*J2 == 0.0){
            continue;
        }
        const double* const J4 = rebx_get_param(rebx, particles[i].ap, "J4");
        const double* const R_eq = rebx_get_param(rebx, particles[i].ap, "R_eq");
        if (R_eq == NULL){
            continue;
        }
        struct reb_vec3d Omega = DEFAULTOMEGA;
        const struct reb_vec3d* Omegaptr = rebx_get_param(rebx, particles[i].ap, "Omega");
        if (Omegaptr != NULL){
            Omega.x = Omegaptr->x;
            Omega.y = Omegaptr->y;
            Omega.z = Omegaptr->z;
        }
        const struct reb_particle pi = particles[i];

        /* new coordinate basis (body-fixed) */
        struct reb_vec3d hatu = {0};
        struct reb_vec3d hatv = {0};
        struct reb_vec3d hatw = {0};

        uvw(Omega, &hatu, &hatv, &hatw);

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
            const double dw = hatw.x*dx + hatw.y*dy + hatw.z*dz;
            const double costheta = dw/r;
            const double costheta2 = costheta*costheta;

            j2_potential_func(G, pi.m, pj.m, J2, R_eq, r, r2, costheta2, &H);
            j4_potential_func(G, pi.m, pj.m, J4, R_eq, r, r2, costheta2, &H);
        }
    }
    return H;
}


