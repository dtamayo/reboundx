/**
 * @file    gravitational_harmonics.c
 * @brief   Add arbitrary-degree gravitational spherical harmonics (J2, J4, ..., tesseral terms) to particles
 * @author  Adaptado de la implementacion original de Dan Tamayo (J2/J4 zonales)
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
 * **Particle Parameters**
 *
 * ============================ ================================= ==================================================================
 * Field (C type)               Required                          Description
 * ============================ ================================= ==================================================================
 * J2 (double)                   No                               J2 coefficient
 * J4 (double)                   No                               J4 coefficient
 * J6 (double)                   No                               J6 coefficient
 * J8 (double)                   No                               J8 coefficient
 * J10 (double)                  No                               J10 coefficient
 * geopotential_model (pointer)  Yes (marks particle as central)  Pointer to rebx_geopotential_model created with rebx_geopotential_create()
 * R_eq (double)                 No                               Equatorial radius of nonspherical body used for calculating Jn harmonics
 * Omega (reb_vec3d)             No                               Angular rotation frequency (Omega_x, Omega_y, Omega_z)
 * theta0 (double)               No                               Length of the body's reference meridian at t=0
 * ============================ ================================= ==================================================================
 *
 * 
 * In difference to the J2/J4 original, this generalized model aims to handle 
 * all the terms needed, included the teseral terms (m > 0). 
 * 
 * =============================================================================
 * HIGHER-ORDER ZONAL HARMONICS (J6, J8, J10)
 * =============================================================================
 * The acceleration vectors (au, av, aw) are obtained by taking the analytical 
 * negative gradient of the zonal geopotential U_n in Cartesian coordinates:
 * 
 *     a_n = -\nabla U_n(x, y, z)
 * 
 * To avoid trigonometric singularities at the poles and eliminate expensive
 * floating-point power operations (pow), the Legendre polynomial derivatives 
 * P'_n(sin(phi)) are expanded explicitly in terms of the normalized coordinate 
 * costheta2 = (z/r)^2. 
 * 
 * The polynomial coefficients (e.g., 429, 12155, 88179) correspond to the exact 
 * algebraic expansions of the Legendre recursive differentials for n = 6, 8, 10.
 * Radial denominator powers are factored into 'f1' to utilize FMA instructions.
 * 
 * Numerical and algebraic integrity verified symbolically via SymPy 3D gradient 
 * differentiation against the canonical Legendre potential definition. Its visible 
 * at the Test_symbolic_zonal_harmonics.py .
 * 
 * 
 * Note: 
 * Effective N is calculated internally after the model is 
 * built. By hard-coding only accepting models with N=50 as maximum. 
 * 
 * References:
 *   - Montenbruck, O., & Gill, E. (2000). Satellite Orbits: Models, Methods and 
 *     Applications. Springer-Verlag. (Sec. 3.2, Eqs. 3.24 - 3.33). 
 * (Not the same model, but we take in consideration the idea of the recursions)
 * ============================================================================= 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"
#include "geopotential.h"

#define DEFAULTOMEGA {0.0, 0.0, 0.0}





static inline void j2_func(double G, double m, const double* J2, const double* R_eq, double r, double r2, double costheta2, double du, double dv, double dw, double* au, double* av, double* aw) {

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

static inline void j4_func(double G, double m, const double* J4, const double* R_eq, double r, double r2, double costheta2, double du, double dv, double dw, double* au, double* av, double* aw) {

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


/**
 * @brief New in this PR: J6 zonal contribution to the acceleration,
 * following the same closed-form pattern as the pre-existing j2_func/
 * j4_func. See the file-level doc comment for the derivation approach
 * (explicit expansion in costheta2 = (z/r)^2, verified symbolically in
 * test_symbolic_geopotential.py). j8_func and j10_func follow the
 * identical pattern for J8 and J10.
 *
 * @param G, m     Standard gravitational parameter pieces (G, central
 *                  body mass).
 * @param J6       Pointer to the particle's "J6" parameter; no-op if
 *                  NULL or zero.
 * @param R_eq     Equatorial radius.
 * @param r, r2    Distance and its square, between bodies.
 * @param costheta2 Square of the cosine of the polar angle (dw/r)^2.
 * @param du, dv, dw  Position offset components in the (u, v, w) frame.
 * @param au, av, aw [in/out] Accumulated acceleration in (u, v, w).
 */

static inline void j6_func(double G, double m, const double* J6, const double* R_eq, double r, double r2, double costheta2, double du, double dv, double dw, double* au, double* av, double* aw) {

    if (J6 == NULL) { return; }
    if (*J6 == 0.0) { return; }

    const double Req2 = (*R_eq)*(*R_eq);
    const double Req6 = Req2*Req2*Req2;
    const double c2 = costheta2*costheta2;
    const double c3 = c2*costheta2;
 
    const double f1 = 7.0/16.0*G*m*(*J6)*Req6/r2/r2/r2/r2/r; 
    const double f2 = 429.0*c3 - 495.0*c2 + 135.0*costheta2 - 5.0;
    const double f3 = 429.0*c3 - 693.0*c2 + 315.0*costheta2 - 35.0;

    *au += f1*f2*du;
    *av += f1*f2*dv;
    *aw += f1*f3*dw;

    return;
}

static inline void j8_func(double G, double m, const double* J8, const double* R_eq, double r, double r2, double costheta2, double du, double dv, double dw, double* au, double* av, double* aw) {

    if (J8 == NULL) { return; }
    if (*J8 == 0.0) { return; }

    const double Req2 = (*R_eq)*(*R_eq);
    const double Req8 = Req2*Req2*Req2*Req2;
    const double c2 = costheta2*costheta2;
    const double c3 = c2*costheta2;
    const double c4 = c2*c2;
 
    const double f1 = 9.0/128.0*G*m*(*J8)*Req8/r2/r2/r2/r2/r2/r; 
    const double f2 = 12155.0*c4 - 20020.0*c3 + 10010.0*c2 - 1540.0*costheta2 + 35.0;
    const double f3 = 12155.0*c4 - 25740.0*c3 + 18018.0*c2 - 4620.0*costheta2 + 315.0;

    *au += f1*f2*du;
    *av += f1*f2*dv;
    *aw += f1*f3*dw;

    return;
}

static inline void j10_func(double G, double m, const double* J10, const double* R_eq, double r, double r2, double costheta2, double du, double dv, double dw, double* au, double* av, double* aw) {

    if (J10 == NULL) { return; }
    if (*J10 == 0.0) { return; }

    const double Req2 = (*R_eq)*(*R_eq);
    const double Req10 = Req2*Req2*Req2*Req2*Req2;
    const double c2 = costheta2*costheta2;
    const double c3 = c2*costheta2;
    const double c4 = c2*c2;
    const double c5 = c4*costheta2;
 
    const double f1 = 11.0/256.0*G*m*(*J10)*Req10/r2/r2/r2/r2/r2/r2/r;
    const double f2 = 88179.0*c5 - 188955.0*c4 + 139230.0*c3 - 40950.0*c2 + 4095.0*costheta2 - 63.0;
    const double f3 = 88179.0*c5 - 230945.0*c4 + 218790.0*c3 - 90090.0*c2 + 15015.0*costheta2 - 693.0;

    *au += f1*f2*du;
    *av += f1*f2*dv;
    *aw += f1*f3*dw;

    return;
}


static inline void j2_potential_func(double G, double mi, double mj, const double* J2, const double* R_eq, double r, double r2, double costheta2, double* H) {

    if (J2 == NULL) { return; }
    if (*J2 == 0.0) { return; }

    const double f1 = G*mi*mj*(*J2)*(*R_eq)*(*R_eq)/r2/r;
    const double f2 = 1.0/2.0*(3.0*costheta2 - 1.0);

    *H += f1*f2;

    return;
}

static inline void j4_potential_func(double G, double mi, double mj, const double* J4, const double* R_eq, double r, double r2, double costheta2, double* H) {

    if (J4 == NULL) { return; }
    if (*J4 == 0.0) { return; }

    const double f1 = G*mi*mj*(*J4)*(*R_eq)*(*R_eq)*(*R_eq)*(*R_eq)/r2/r2/r;
    const double f2 = 1.0/8.0*(35.0*costheta2*costheta2 - 30.0*costheta2 + 3.0);

    *H += f1*f2;

    return;
}


static inline void j6_potential_func(double G, double mi, double mj, const double* J6, const double* R_eq, double r, double r2, double costheta2, double* H) {
 
    if (J6 == NULL) { return; }
    if (*J6 == 0.0) { return; }
 
    const double Req2 = (*R_eq)*(*R_eq);
    const double Req6 = Req2*Req2*Req2;
    const double c2 = costheta2*costheta2;
    const double c3 = c2*costheta2;
 
    const double f1 = G*mi*mj*(*J6)*Req6/r2/r2/r2/r; 
    const double f2 = 1.0/16.0*(231.0*c3 - 315.0*c2 + 105.0*costheta2 - 5.0);
 
    *H += f1*f2;
 
    return;
}
 
static inline void j8_potential_func(double G, double mi, double mj, const double* J8, const double* R_eq, double r, double r2, double costheta2, double* H) {
 
    if (J8 == NULL) { return; }
    if (*J8 == 0.0) { return; }
 
    const double Req2 = (*R_eq)*(*R_eq);
    const double Req8 = Req2*Req2*Req2*Req2;
    const double c2 = costheta2*costheta2;
    const double c3 = c2*costheta2;
    const double c4 = c2*c2;
 
    const double f1 = G*mi*mj*(*J8)*Req8/r2/r2/r2/r2/r; 
    const double f2 = 1.0/128.0*(6435.0*c4 - 12012.0*c3 + 6930.0*c2 - 1260.0*costheta2 + 35.0);
 
    *H += f1*f2;
 
    return;
}
 
static inline void j10_potential_func(double G, double mi, double mj, const double* J10, const double* R_eq, double r, double r2, double costheta2, double* H) {
 
    if (J10 == NULL) { return; }
    if (*J10 == 0.0) { return; }
 
    const double Req2 = (*R_eq)*(*R_eq);
    const double Req10 = Req2*Req2*Req2*Req2*Req2;
    const double c2 = costheta2*costheta2;
    const double c3 = c2*costheta2;
    const double c4 = c2*c2;
    const double c5 = c4*costheta2;
 
    const double f1 = G*mi*mj*(*J10)*Req10/r2/r2/r2/r2/r2/r; 
    const double f2 = 1.0/256.0*(46189.0*c5 - 109395.0*c4 + 90090.0*c3 - 30030.0*c2 + 3465.0*costheta2 - 63.0);
 
    *H += f1*f2;
 
    return;
}

/**
 * @brief Build an orthonormal frame (hatu, hatv, hatw) with hatw aligned
 * with the rotation axis Omega (or the z-axis if Omega is zero).
 *
 * @param Omega Angular velocity vector of the central body.
 * @param hatu, hatv, hatw [out] Orthonormal right-handed frame.
 */

static inline void uvw(struct reb_vec3d Omega, struct reb_vec3d* hatu, struct reb_vec3d* hatv, struct reb_vec3d* hatw) {

    const double omega2 = Omega.x*Omega.x + Omega.y*Omega.y + Omega.z*Omega.z;
    const double omega = sqrt(omega2);

    struct reb_vec3d s = {0, 0, 1.0};
    if (omega > 0.0) {
        s.x = Omega.x/omega;
        s.y = Omega.y/omega;
        s.z = Omega.z/omega;
    }

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

/**
 * @brief Read a particle's "Omega" (angular velocity vector) and
 * "theta0" (reference meridian angle at t=0) parameters, defaulting to
 * a non-rotating body if absent, and return the current rotation angle
 * and the (u, v, w) frame aligned with the rotation axis.
 *
 * @param pi   Central body particle.
 * @param t    Current simulation time.
 * @param hatu, hatv, hatw [out] Orthonormal frame from uvw(Omega, ...).
 * @return Current rotation angle theta0 + |Omega| * t, radians.
 */

static double rebx_geopotential_read_orientation(struct rebx_extras* const rebx, const struct reb_particle* pi,
                                                  double t, struct reb_vec3d *hatu, struct reb_vec3d *hatv, struct reb_vec3d *hatw) {
    struct reb_vec3d Omega = DEFAULTOMEGA;
    const struct reb_vec3d* Omegaptr = rebx_get_param(rebx, pi->ap, "Omega");
    if (Omegaptr != NULL){
        Omega = *Omegaptr;
    }
    const double* theta0ptr = rebx_get_param(rebx, pi->ap, "theta0");
    const double theta0 = (theta0ptr != NULL) ? *theta0ptr : 0.0;
 
    const double omega_mag = sqrt(Omega.x*Omega.x + Omega.y*Omega.y + Omega.z*Omega.z);
 
    uvw(Omega, hatu, hatv, hatw);
    return theta0 + omega_mag * t;
}

/**
 * @brief Body-fixed spherical coordinates (phi, r, lambda_body) of pj as
 * seen from pi, expressed in pi's rotating (u, v, w) frame.
 *
 * @param pi, pj   Central body and orbiting particle.
 * @param hatu, hatv, hatw  Orthonormal frame aligned with pi's rotation
 *                          axis (hatw), from uvw().
 * @param theta    Current rotation angle of the body's reference
 *                 meridian about hatw (theta0 + |Omega| * t), from
 *                 rebx_geopotential_read_orientation().
 * @param phi         [out] Body-fixed latitude, radians.
 * @param r           [out] Distance between pi and pj.
 * @param lambda_body [out] Body-fixed (co-rotating) longitude, radians.
 */

static void rebx_geopotential_geometry(const struct reb_particle* pi, const struct reb_particle* pj,
                                        struct reb_vec3d hatu, struct reb_vec3d hatv, struct reb_vec3d hatw,
                                        double theta,
                                        double *phi, double *r, double *lambda_body) {
    const double dx = pj->x - pi->x;
    const double dy = pj->y - pi->y;
    const double dz = pj->z - pi->z;
 
    const double du = hatu.x*dx + hatu.y*dy + hatu.z*dz;
    const double dv = hatv.x*dx + hatv.y*dy + hatv.z*dz;
    const double dw = hatw.x*dx + hatw.y*dy + hatw.z*dz;
 
    *r = sqrt(du*du + dv*dv + dw*dw);
    *phi = asin(dw / (*r));
    *lambda_body = atan2(dv, du) - theta;
}


/**
 * @brief Convert (a_r, a_phi, a_lambda) spherical acceleration components
 * at a given (sinphi, cosphi, sinlam, coslam) into Cartesian components
 * along the (u, v, w) frame.
 *
 * @param a_r, a_phi, a_lambda  Spherical acceleration components, from
 *                              rebx_geopotential_acceleration().
 * @param sinphi, cosphi  sin/cos of the body-fixed latitude.
 * @param sinlam, coslam  sin/cos of the longitude in the (u, v) plane
 *                        (inertial-frame longitude, i.e. lambda_body +
 *                        theta — NOT the body-fixed lambda_body used to
 *                        evaluate the model).
 * @param au, av, aw [out] Acceleration components along (hatu, hatv, hatw).
 */

static inline void spherical_to_uvw(double a_r, double a_phi, double a_lambda,
                                     double sinphi, double cosphi,
                                     double sinlam, double coslam,
                                     double *au, double *av, double *aw) {


    const double er[3]      = { cosphi*coslam,  cosphi*sinlam,  sinphi };
    const double ephi[3]    = { -sinphi*coslam, -sinphi*sinlam, cosphi };
    const double elambda[3] = { -sinlam,          coslam,         0.0  };

    *au = a_r*er[0] + a_phi*ephi[0] + a_lambda*elambda[0];
    *av = a_r*er[1] + a_phi*ephi[1] + a_lambda*elambda[1];
    *aw = a_r*er[2] + a_phi*ephi[2] + a_lambda*elambda[2];
}


void rebx_gravitational_harmonics(struct reb_simulation* const sim, struct rebx_force* const gh, struct reb_particle* const particles, const int N){
    const double G = sim->G;
    const double t = sim->t;
    struct rebx_extras* const rebx = sim->extras;

    for (int i=0; i<N; i++){

        rebx_geopotential_model* const model = (rebx_geopotential_model*) rebx_get_param(rebx, particles[i].ap, "geopotential_model");
        if (model != NULL) {
            /* ---------- general root---------- */

            const struct reb_particle pi = particles[i];
 
            struct reb_vec3d hatu = {0}, hatv = {0}, hatw = {0};
            const double theta = rebx_geopotential_read_orientation(rebx, &pi, t, &hatu, &hatv, &hatw);
 
            const double costheta = cos(theta);
            const double sintheta = sin(theta);

            for (int j=0; j<N; j++){
                if (j == i){
                    continue;
                }
                const struct reb_particle pj = particles[j];
 
                double phi, r, lambda_body;
                rebx_geopotential_geometry(&pi, &pj, hatu, hatv, hatw, theta, &phi, &r, &lambda_body);
 
                double a_r, a_phi, a_lambda;
                double sinphi, cosphi;
                rebx_geopotential_acceleration(model, phi, r, lambda_body, &a_r, &a_phi, &a_lambda, &sinphi, &cosphi);
 
                const double coslam_body = model->cosml[1];
                const double sinlam_body = model->sinml[1];
                const double coslam = coslam_body * costheta - sinlam_body * sintheta;
                const double sinlam = sinlam_body * costheta + coslam_body * sintheta;

                double au, av, aw;
                spherical_to_uvw(a_r, a_phi, a_lambda, sinphi, cosphi, sinlam, coslam, &au, &av, &aw);
 
                const double ax = hatu.x*au + hatv.x*av + hatw.x*aw;
                const double ay = hatu.y*au + hatv.y*av + hatw.y*aw;
                const double az = hatu.z*au + hatv.z*av + hatw.z*aw;
 
                particles[j].ax += ax;
                particles[j].ay += ay;
                particles[j].az += az;
 
                const double fac = pj.m/pi.m;
                particles[i].ax -= fac*ax;
                particles[i].ay -= fac*ay;
                particles[i].az -= fac*az;
            }
            continue;   
        }
        /* ---------- previous root (J2/J4, retrocompatible) ---------- */

        const double* const J2 = rebx_get_param(rebx, particles[i].ap, "J2");
        if (J2 == NULL){
            continue;
        }
        if (*J2 == 0.0){
            continue;
        }
        const double* const J4 = rebx_get_param(rebx, particles[i].ap, "J4");
        const double* const J6 = rebx_get_param(rebx, particles[i].ap, "J6");
        const double* const J8 = rebx_get_param(rebx, particles[i].ap, "J8");
        const double* const J10 = rebx_get_param(rebx, particles[i].ap, "J10");
        const double* const R_eq = rebx_get_param(rebx, particles[i].ap, "R_eq");
        if (R_eq == NULL){
            continue;
        }

        struct reb_vec3d Omega = DEFAULTOMEGA;
        const struct reb_vec3d* Omegaptr = rebx_get_param(rebx, particles[i].ap, "Omega");
        if (Omegaptr != NULL){
            Omega = *Omegaptr;
        }

        const struct reb_particle pi = particles[i];

        struct reb_vec3d hatu = {0}, hatv = {0}, hatw = {0};
        uvw(Omega, &hatu, &hatv, &hatw);


        struct reb_vec3d hatx_ = { hatu.x, hatv.x, hatw.x };
        struct reb_vec3d haty_ = { hatu.y, hatv.y, hatw.y };
        struct reb_vec3d hatz_ = { hatu.z, hatv.z, hatw.z };

        for (int j=0; j<N; j++){
            if (j == i){
                continue;
            }
            const struct reb_particle pj = particles[j];
 
            const double dx = pj.x - pi.x;
            const double dy = pj.y - pi.y;
            const double dz = pj.z - pi.z;
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double r  = sqrt(r2);
 
            const double du = hatu.x*dx + hatu.y*dy + hatu.z*dz;
            const double dv = hatv.x*dx + hatv.y*dy + hatv.z*dz;
            const double dw = hatw.x*dx + hatw.y*dy + hatw.z*dz;
            const double costheta  = dw/r;
            const double costheta2 = costheta*costheta;
 
            double au = 0.0, av = 0.0, aw = 0.0;
 
            j2_func(G, pi.m, J2, R_eq, r, r2, costheta2, du, dv, dw, &au, &av, &aw);
            j4_func(G, pi.m, J4, R_eq, r, r2, costheta2, du, dv, dw, &au, &av, &aw);
            j6_func(G, pi.m, J6, R_eq, r, r2, costheta2, du, dv, dw, &au, &av, &aw);
            j8_func(G, pi.m, J8, R_eq, r, r2, costheta2, du, dv, dw, &au, &av, &aw);
            j10_func(G, pi.m, J10, R_eq, r, r2, costheta2, du, dv, dw, &au, &av, &aw);
 
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

double rebx_gravitational_harmonics_potential(struct rebx_extras* const rebx){
    if (rebx->sim == NULL){
        rebx_error(rebx, "");
        return 0;
    }
    struct reb_simulation* const sim = rebx->sim;
    struct reb_particle* const particles = sim->particles;
    const double G = sim->G;
    const int N = sim->N - sim->N_var;
    double H = 0.0;
 
    for (int i=0; i<N; i++){
 
        rebx_geopotential_model* const model = rebx_get_param(rebx, particles[i].ap, "geopotential_model");
 
        if (model != NULL){
            /* ---------- ruta general ---------- */
            
            const struct reb_particle pi = particles[i];
 
            struct reb_vec3d hatu = {0}, hatv = {0}, hatw = {0};
            const double theta = rebx_geopotential_read_orientation(rebx, &pi, sim->t, &hatu, &hatv, &hatw);
 
            for (int j=0; j<N; j++){
                if (j == i){
                    continue;
                }
                const struct reb_particle pj = particles[j];
 
                double phi, r, lambda_body;
                rebx_geopotential_geometry(&pi, &pj, hatu, hatv, hatw, theta, &phi, &r, &lambda_body);
 
                const double U_pert = rebx_geopotential_potential(model, phi, r, lambda_body);
                H -= pj.m * U_pert;
            }
            continue;
        }
 
        /* ---------- ruta legada (J2/J4) ---------- */
        const double* const J2 = rebx_get_param(rebx, particles[i].ap, "J2");
        if (J2 == NULL){
            continue;
        }
        if (*J2 == 0.0){
            continue;
        }
        
        const double* const J4 = rebx_get_param(rebx, particles[i].ap, "J4");
        const double* const J6 = rebx_get_param(rebx, particles[i].ap, "J6");
        const double* const J8 = rebx_get_param(rebx, particles[i].ap, "J8");
        const double* const J10 = rebx_get_param(rebx, particles[i].ap, "J10");
        const double* const R_eq = rebx_get_param(rebx, particles[i].ap, "R_eq");

        if (R_eq == NULL){
            continue;
        }
        struct reb_vec3d Omega = DEFAULTOMEGA;
        const struct reb_vec3d* Omegaptr = rebx_get_param(rebx, particles[i].ap, "Omega");
        if (Omegaptr != NULL){
            Omega = *Omegaptr;
        }
        const struct reb_particle pi = particles[i];
 
        struct reb_vec3d hatu = {0}, hatv = {0}, hatw = {0};
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
            const double r  = sqrt(r2);
 
            const double dw = hatw.x*dx + hatw.y*dy + hatw.z*dz;
            const double costheta  = dw/r;
            const double costheta2 = costheta*costheta;
 
            j2_potential_func(G, pi.m, pj.m, J2, R_eq, r, r2, costheta2, &H);
            j4_potential_func(G, pi.m, pj.m, J4, R_eq, r, r2, costheta2, &H);
            j6_potential_func(G, pi.m, pj.m, J6, R_eq, r, r2, costheta2, &H);
            j8_potential_func(G, pi.m, pj.m, J8, R_eq, r, r2, costheta2, &H);
            j10_potential_func(G, pi.m, pj.m, J10, R_eq, r, r2, costheta2, &H);

        }
    }
    return H;
}