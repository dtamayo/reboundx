/**
 * @file    triaxial_torque.c
 * @brief   Torque on triaxial bodies
 * @author  Henry Yuan
 *
 *
 * *** COMMENT SECTION BELOW THIS HAS NOT BEEN CHANGED FROM COPIED GRAVITATIONAL_HARMONICS.C
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
 * Authors                 D. Tamayo
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
 * Based on                None
 * C Example               :ref:`c_example_J2`
 * Python Example          `J2.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/J2.ipynb>`_.
 * ======================= ===============================================
 *
 * * *** COMMENT SECTION ABOVE THIS HAS NOT BEEN CHANGED FROM COPIED GRAVITATIONAL_HARMONICS.C
 *
 * Adds the effects of a particle having 3 differing moments of inertia (triaxial particle) on the spin vector
 * of the particle. Assumes that changes to the spin angular momentum are negligible compared to the orbital
 * angular momentum; i.e., the particle's orbit is not affected by the changes to the spin vector.
 *
 *
 * **Effect Parameters**
 *
 * None
 *
 * **Particle Parameters**
 *
 * x,y,z: refers to the x,y,z coordinate system underlying rebound
 *
 * I define a new coordinate system in addition to the one used by rebound (described above):
 * -> i,j,k: principal axes corresponding to the particle's principal moments of inertia. axis i has the
 *    lowest moment of inertia and axis k has the highest moment of inertia
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * tt_Ii (double)                  Yes         Moment for axis i (<= Ij)
 * tt_Ij (double)                  Yes         Moment for axis j (<= Ik, >= Ii)
 * tt_Ik (double)                  Yes         Moment for axis k (>= Ij)
 * tt_omega (double)               Yes         spin rate
 * tt_ix (double)                  Yes         x-component of unit vector for lowest moment axis
 * tt_iy (double)                  Yes         y-component of unit vector for lowest moment axis
 * tt_iz (double)                  Yes         z-component of unit vector for lowest moment axis
 * tt_jx (double)                  Yes         x-component of unit vector for middle moment axis
 * tt_jy (double)                  Yes         y-component of unit vector for middle moment axis
 * tt_jz (double)                  Yes         z-component of unit vector for middle moment axis
 * tt_kx (double)                  Yes         x-component of unit vector for highest moment axis
 * tt_ky (double)                  Yes         y-component of unit vector for highest moment axis
 * tt_kz (double)                  Yes         z-component of unit vector for highest moment axis
 * tt_si (double)                  Yes         i component of spin unit vector
 * tt_sj (double)                  Yes         j component of spin unit vector
 * tt_sk (double)                  Yes         k component of spin unit vector
 * ============================ =========== ==================================================================
 *
 * Parameter reuqirements:
 * ============================
 * i, j, k, and s must be unit vectors (e.g. such that ix^2 + iy^2 + iz^2 = 1)
 *
 */

/* DEBUGGING CHANGES:
- commented out torque calc
- domega_dts all 0
- print statements
- only using first derivative calc
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

// linearly interpolates position of particle p at time dt
static void rebx_interpolate_xyz(struct reb_particle* p, double xyz[3], const double dt){
    xyz[0] = p->x + (dt*p->vx);
    xyz[1] = p->y + (dt*p->vy);
    xyz[2] = p->z + (dt*p->vz);
}

// convert vector ijk to xyz (same vector in xyz basis), given ijk_xyz
static void rebx_ijk_to_xyz(double ijk[3], double xyz[3], double ijk_xyz[3][3]){
    for (int i=0; i<3; i++) {
        xyz[i] = ijk[0]*ijk_xyz[0][i] + ijk[1]*ijk_xyz[1][i] + ijk[2]*ijk_xyz[2][i];
    }
}

// computes time-derivative of spin vector omega_i,omega_j,omega_k using torque vector Mi, Mj, Mk according to Euler's equations
static void rebx_domega_dt(double omega_ijk[3], double M_ijk[3], double I_ijk[3], double domega_dts[3]){
    domega_dts[0] = (M_ijk[0] + (I_ijk[1]-I_ijk[2])*omega_ijk[1]*omega_ijk[2]) / I_ijk[0];
    domega_dts[1] = (M_ijk[1] + (I_ijk[2]-I_ijk[0])*omega_ijk[2]*omega_ijk[0]) / I_ijk[1];
    domega_dts[2] = (M_ijk[2] + (I_ijk[0]-I_ijk[1])*omega_ijk[0]*omega_ijk[1]) / I_ijk[2];
    // domega_dts[0] = 0.0; // [DEBUG]
    // domega_dts[1] = 0.0; // [DEBUG]
    // domega_dts[2] = 0.0; // [DEBUG]
}

// computes time-derivative of vectors i,j,k (components in old ijk basis)
static void rebx_dijk_dt(double ijk_ijk[3][3], double omega_ijk[3], double dijk_dts[3][3]){
    // di/dt
    dijk_dts[0][0] = omega_ijk[1]*ijk_ijk[0][2] - omega_ijk[2]*ijk_ijk[0][1];
    dijk_dts[0][1] = omega_ijk[2]*ijk_ijk[0][0] - omega_ijk[0]*ijk_ijk[0][2];
    dijk_dts[0][2] = omega_ijk[0]*ijk_ijk[0][1] - omega_ijk[1]*ijk_ijk[0][0];

    // dj/dt
    dijk_dts[1][0] = omega_ijk[1]*ijk_ijk[1][2] - omega_ijk[2]*ijk_ijk[1][1];
    dijk_dts[1][1] = omega_ijk[2]*ijk_ijk[1][0] - omega_ijk[0]*ijk_ijk[1][2];
    dijk_dts[1][2] = omega_ijk[0]*ijk_ijk[1][1] - omega_ijk[1]*ijk_ijk[1][0];

    // dk/dt
    dijk_dts[2][0] = omega_ijk[1]*ijk_ijk[2][2] - omega_ijk[2]*ijk_ijk[2][1];
    dijk_dts[2][1] = omega_ijk[2]*ijk_ijk[2][0] - omega_ijk[0]*ijk_ijk[2][2];
    dijk_dts[2][2] = omega_ijk[0]*ijk_ijk[2][1] - omega_ijk[1]*ijk_ijk[2][0];
}

// calculates the triaxial torque from all other bodies on the 'index'th particle
static void rebx_calc_torques(struct reb_simulation* const sim, int index, double M_ijk[3], double I_ijk[3], double ijk_xyz[3][3], const double dt, const double dt_tot){

    struct reb_particle* p = &sim->particles[index];
    struct reb_particle* torquer;
    double p_xyz[3];
    double torquer_xyz[3];
    double rx;
    double ry;
    double rz;
    double r;
    double r_dot_i;
    double r_dot_j;
    double r_dot_k;
    double prefac;

    rebx_interpolate_xyz(p,p_xyz,dt - dt_tot);
    M_ijk[0] = 0;
    M_ijk[1] = 0;
    M_ijk[2] = 0;

    const int _N_real = sim->N - sim->N_var;
    for(int i=0; i<_N_real; i++){
        if (i == index) {
            continue;
        }
        torquer = &sim->particles[i];
        rebx_interpolate_xyz(torquer,torquer_xyz,dt - dt_tot);
        rx = p_xyz[0] - torquer_xyz[0];
        ry = p_xyz[1] - torquer_xyz[1];
        rz = p_xyz[2] - torquer_xyz[2];
        r = sqrt(rx*rx + ry*ry + rz*rz);
        prefac = 3 * sim->G * torquer->m / pow(r,5);

        r_dot_i = rx*ijk_xyz[0][0] + ry*ijk_xyz[0][1] + rz*ijk_xyz[0][2];
        r_dot_j = rx*ijk_xyz[1][0] + ry*ijk_xyz[1][1] + rz*ijk_xyz[1][2];
        r_dot_k = rx*ijk_xyz[2][0] + ry*ijk_xyz[2][1] + rz*ijk_xyz[2][2];

        M_ijk[0] += prefac*(I_ijk[2]-I_ijk[1])*r_dot_j*r_dot_k;
        M_ijk[1] += prefac*(I_ijk[0]-I_ijk[2])*r_dot_k*r_dot_i;
        M_ijk[2] += prefac*(I_ijk[1]-I_ijk[0])*r_dot_i*r_dot_j;
    }
}

/* updates spin vector, omega, and ijk in lockstep using 4th order Runge Kutta.
If calc_torque_bool = 0, torque NOT calculated, otherwise torque calculated */
static void rebx_update_spin_ijk(struct reb_simulation* const sim, int calc_torque_bool, int index, double* const ix, double* const iy, double* const iz,
    double* const jx, double* const jy, double* const jz, double* const kx, double* const ky, double* const kz, double* const si,
    double* const sj, double* const sk, double* const omega, const double Ii, const double Ij, const double Ik, const double dt){

    // Array for principal moments
    double I_ijk[3];
    I_ijk[0] = Ii;
    I_ijk[1] = Ij;
    I_ijk[2] = Ik;

    // Declare matrices for all R-K calculations
    double rk_M_ijk[4][3] = {}; // ijk components of torque on body, Needs all values initialized to zero because multiple functions add to the values
    double rk_omega_ijk[4][3]; // ijk components of spin vector (total omega vector)
    double rk_ijk_xyz[4][3][3]; // xyz components of each ijk vector
    double rk_ijk_ijk[4][3][3]; // components of each ijk vector in the old ijk basis
    double rk_domega_dts[4][3]; // matrix for all calculations of domega/dt
    double rk_dijk_dts[4][3][3]; // matrix for all calculations of d{ijk}/dt
    /* second dimension: vector (i_hat,j_hat,k_hat)
    third dimension: component of vector (i,j,k) */
    double rk_dts[4] = {0.0, 0.5*dt, 0.5*dt, dt}; // array of sub-timesteps for each RK calculation

    // initialize omega_0
    rk_omega_ijk[0][0] = *omega**si;
    rk_omega_ijk[0][1] = *omega**sj;
    rk_omega_ijk[0][2] = *omega**sk;

    // initialize xyz components of i,j,k
    rk_ijk_xyz[0][0][0] = *ix;
    rk_ijk_xyz[0][0][1] = *iy;
    rk_ijk_xyz[0][0][2] = *iz;
    rk_ijk_xyz[0][1][0] = *jx;
    rk_ijk_xyz[0][1][1] = *jy;
    rk_ijk_xyz[0][1][2] = *jz;
    rk_ijk_xyz[0][2][0] = *kx;
    rk_ijk_xyz[0][2][1] = *ky;
    rk_ijk_xyz[0][2][2] = *kz;

    // initialize ijk components of i,j,k
    rk_ijk_ijk[0][0][0] = 1.0;
    rk_ijk_ijk[0][0][1] = 0.0;
    rk_ijk_ijk[0][0][2] = 0.0;
    rk_ijk_ijk[0][1][0] = 0.0;
    rk_ijk_ijk[0][1][1] = 1.0;
    rk_ijk_ijk[0][1][2] = 0.0;
    rk_ijk_ijk[0][2][0] = 0.0;
    rk_ijk_ijk[0][2][1] = 0.0;
    rk_ijk_ijk[0][2][2] = 1.0;

    /****************************************************************************************/
    // Runge-Kutta Calculations
    for (int i = 0; i < 4; i++) {
        // Pre-calcs, skip if first iteration
        if (i != 0) {
            for (int j=0; j < 3; j++) {
                rk_omega_ijk[i][j] = rk_domega_dts[i-1][j]*rk_dts[i] + rk_omega_ijk[0][j];
                for (int k=0; k < 3; k++) {
                    rk_ijk_ijk[i][j][k] = rk_dijk_dts[i-1][j][k]*rk_dts[i] + rk_ijk_ijk[0][j][k];
                }
                rebx_ijk_to_xyz(rk_ijk_ijk[i][j],rk_ijk_xyz[i][j],rk_ijk_xyz[0]);
            }
        }

        // Calcs
        if (calc_torque_bool != 0) {
            rebx_calc_torques(sim,index,rk_M_ijk[i],I_ijk,rk_ijk_xyz[i],rk_dts[i], dt); // [DEBUG]
        }
        rebx_domega_dt(rk_omega_ijk[i],rk_M_ijk[i],I_ijk,rk_domega_dts[i]);
        rebx_dijk_dt(rk_ijk_ijk[i],rk_omega_ijk[i],rk_dijk_dts[i]);
    }

    /****************************************************************************************/
    /****************************************************************************************/
    // // First Runge-Kutta calculations
    // rebx_calc_torques(sim,index,M_ijk[0],I_ijk,rk_ijk_xyz[0],0.0);
    // rebx_domega_dt(rk_omega_ijk[0],rk_M_ijk[0],I_ijk,rk_domega_dts[0]);
    // rebx_dijk_dt(rk_ijk_ijk[0],rk_omega_ijk[0],rk_dijk_dts[0]);

    // // Pre-second RK calcs
    // for (int j=0; j < 3; j++) {
    //     rk_omega_ijk[1][j] = rk_domega_dts[0][j]*dt*0.5 + rk_omega_ijk[0][j];
    //     for (int k=0; k < 3; k++) {
    //         rk_ijk_ijk[1][j][k] = rk_dijk_dts[0][j][k]*dt*0.5 + rk_ijk_ijk[0][j][k];
    //     }
    //     rebx_ijk_to_xyz(rk_ijk_ijk[1][j],rk_ijk_xyz[1][j],rk_ijk_xyz[0]);
    // }

    // // Second RK calcs
    // rebx_calc_torques(sim,index,M_ijk[1],I_ijk,rk_ijk_xyz[1],0.5*dt);
    // rebx_domega_dt(rk_omega_ijk[1],rk_M_ijk[1],I_ijk,rk_domega_dts[1]);
    // rebx_dijk_dt(rk_ijk_ijk[1],rk_omega_ijk[1],rk_dijk_dts[1]);

    // Continue editing here
    /****************************************************************************************/
    /****************************************************************************************/

    // calculate domega, d{ijk}
    double domega[3];
    double dijk[3][3];
    for (int i = 0; i < 3; i++){
        // domega[i] = rk_domega_dts[0][i] * dt; // [DEBUG]
        domega[i] = (rk_domega_dts[0][i] + 2*rk_domega_dts[1][i] + 2*rk_domega_dts[2][i] + rk_domega_dts[3][i]) * dt / 6;
        for (int j = 0; j < 3; j++){
            // dijk[i][j] = rk_dijk_dts[0][i][j] * dt; // [DEBUG]
            dijk[i][j] = (rk_dijk_dts[0][i][j] + 2*rk_dijk_dts[1][i][j] + 2*rk_dijk_dts[2][i][j] + rk_dijk_dts[3][i][j]) * dt / 6;
        }
    }

    // printf("%f\n", domega[0]); // [DEBUG]
    // printf("%f\n", domega[1]); // [DEBUG]
    // printf("%f\n", domega[2]); // [DEBUG]

    double omega_i = rk_omega_ijk[0][0] + domega[0];
    double omega_j = rk_omega_ijk[0][1] + domega[1];
    double omega_k = rk_omega_ijk[0][2] + domega[2];

    *omega = sqrt(omega_i*omega_i + omega_j*omega_j + omega_k*omega_k);
    *si = omega_i / *omega;
    *sj = omega_j / *omega;
    *sk = omega_k / *omega;

    double ijk_ijk[3][3];
    double ijk_xyz[3][3];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ijk_ijk[i][j] = rk_ijk_ijk[0][i][j] + dijk[i][j];
        }
        rebx_ijk_to_xyz(ijk_ijk[i],ijk_xyz[i],rk_ijk_xyz[0]);
    }

    double ix_ = ijk_xyz[0][0];
    double iy_ = ijk_xyz[0][1];
    double iz_ = ijk_xyz[0][2];
    double jx_ = ijk_xyz[1][0];
    double jy_ = ijk_xyz[1][1];
    double jz_ = ijk_xyz[1][2];
    double kx_ = ijk_xyz[2][0];
    double ky_ = ijk_xyz[2][1];
    double kz_ = ijk_xyz[2][2];

    // re-normalize
    double i_mag = sqrt(ix_*ix_ + iy_*iy_ + iz_*iz_);
    double j_mag = sqrt(jx_*jx_ + jy_*jy_ + jz_*jz_);
    double k_mag = sqrt(kx_*kx_ + ky_*ky_ + kz_*kz_);

    *ix = ix_ / i_mag;
    *iy = iy_ / i_mag;
    *iz = iz_ / i_mag;
    *jx = jx_ / j_mag;
    *jy = jy_ / j_mag;
    *jz = jz_ / j_mag;
    *kx = kx_ / k_mag;
    *ky = ky_ / k_mag;
    *kz = kz_ / k_mag;
}

/* START YS CHANGES */
static double mag_vec(const double vec[3]) {
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}
static void renorm_vec(double vec[3]) {
    double vec_mag = mag_vec(vec);
    vec[0] /= vec_mag;
    vec[1] /= vec_mag;
    vec[2] /= vec_mag;
}
static void renorm_mat(double vec[3][3]) {
    renorm_vec(vec[0]);
    renorm_vec(vec[1]);
    renorm_vec(vec[2]);
}
static void rebx_update_spin_ijk_euler(struct reb_simulation* const sim, int index, double* const ix, double* const iy, double* const iz,
    double* const jx, double* const jy, double* const jz, double* const kx, double* const ky, double* const kz, double* const si,
    double* const sj, double* const sk, double* const omega, const double Ii, const double Ij, const double Ik, const double dt){
    // Array for principal moments
    double I_ijk[3];
    I_ijk[0] = Ii;
    I_ijk[1] = Ij;
    I_ijk[2] = Ik;

    double M_ijk[3]; // ijk components of torque on body. Let "calc torques" set to zero
    double ijk_xyz[3][3]; // xyz components of each ijk vector
    double ijk_ijk[3][3]; // components of each ijk vector in the old ijk basis
    double domega_dts[3]; // matrix for all calculations of domega/dt
    double dijk_dts[3][3]; // matrix for all calculations of d{ijk}/dt
    double omega_ijk[3];

    /****** INITIALIZATIONS ******/
    ijk_xyz[0][0] = *ix;
    ijk_xyz[0][1] = *iy;
    ijk_xyz[0][2] = *iz;
    ijk_xyz[1][0] = *jx;
    ijk_xyz[1][1] = *jy;
    ijk_xyz[1][2] = *jz;
    ijk_xyz[2][0] = *kx;
    ijk_xyz[2][1] = *ky;
    ijk_xyz[2][2] = *kz;
    ijk_ijk[0][0] = 1.0;
    ijk_ijk[0][1] = 0.0;
    ijk_ijk[0][2] = 0.0;
    ijk_ijk[1][0] = 0.0;
    ijk_ijk[1][1] = 1.0;
    ijk_ijk[1][2] = 0.0;
    ijk_ijk[2][0] = 0.0;
    ijk_ijk[2][1] = 0.0;
    ijk_ijk[2][2] = 1.0;
    omega_ijk[0] = *omega**si;
    omega_ijk[1] = *omega**sj;
    omega_ijk[2] = *omega**sk;

    /******* CALCULATE all dY/dt (Y = each variable) ***************************/
    rebx_calc_torques(sim,index,M_ijk,I_ijk,ijk_xyz,0, dt);
    rebx_domega_dt(omega_ijk,M_ijk,I_ijk,domega_dts);
    rebx_dijk_dt(ijk_ijk,omega_ijk,dijk_dts);

    /******* CALCULATE all dYs ***************************/
    double domega[3];
    double dijk[3][3];
    for (int i = 0; i < 3; i++){
        domega[i] = domega_dts[i] * dt;
        for (int j = 0; j < 3; j++){
            dijk[i][j] = dijk_dts[i][j] * dt;
        }
    }

    /******* CALCULATE all new Y's ***************************/
    double new_omega[3];
    double new_ijk_ijk[3][3];
    double new_ijk_xyz[3][3];

    for (int i = 0; i < 3; i++) {
        new_omega[i] = omega_ijk[i] + domega[i];
        for (int j = 0; j < 3; j++) {
            new_ijk_ijk[i][j] = ijk_ijk[i][j] + dijk[i][j];
        }
        rebx_ijk_to_xyz(new_ijk_ijk[i],new_ijk_xyz[i],ijk_xyz);
    }
    renorm_mat(new_ijk_xyz);

    /******* UPDATE ALL Y's ***************************/
    *omega = mag_vec(new_omega);
    renorm_vec(new_omega);
    *si = new_omega[0];
    *sj = new_omega[1];
    *sk = new_omega[2];
    *ix = new_ijk_xyz[0][0];
    *iy = new_ijk_xyz[0][1];
    *iz = new_ijk_xyz[0][2];
    *jx = new_ijk_xyz[1][0];
    *jy = new_ijk_xyz[1][1];
    *jz = new_ijk_xyz[1][2];
    *kx = new_ijk_xyz[2][0];
    *ky = new_ijk_xyz[2][1];
    *kz = new_ijk_xyz[2][2];
}

static void rebx_update_spin_ijk_mid(struct reb_simulation* const sim, int index, double* const ix, double* const iy, double* const iz,
    double* const jx, double* const jy, double* const jz, double* const kx, double* const ky, double* const kz, double* const si,
    double* const sj, double* const sk, double* const omega, const double Ii, const double Ij, const double Ik, const double dt){
    // Array for principal moments
    double I_ijk[3];
    I_ijk[0] = Ii;
    I_ijk[1] = Ij;
    I_ijk[2] = Ik;

    double M_ijk[3]; // ijk components of torque on body. Let "calc torques" set to zero
    double ijk_xyz[3][3]; // xyz components of each ijk vector
    double ijk_ijk[3][3]; // components of each ijk vector in the old ijk basis
    double domega_dts[3]; // matrix for all calculations of domega/dt
    double dijk_dts[3][3]; // matrix for all calculations of d{ijk}/dt
    double omega_ijk[3];

    ijk_xyz[0][0] = *ix;
    ijk_xyz[0][1] = *iy;
    ijk_xyz[0][2] = *iz;
    ijk_xyz[1][0] = *jx;
    ijk_xyz[1][1] = *jy;
    ijk_xyz[1][2] = *jz;
    ijk_xyz[2][0] = *kx;
    ijk_xyz[2][1] = *ky;
    ijk_xyz[2][2] = *kz;

    /******** INITIALIZE Y(t), dYdt(t), dY FOR EACH VAR (2x) ******/
    double domega[3];
    double dijk[3][3];
    /* we have a small trick here: for a N-step method, store (N+1) copies of
     * each dynamical variable.
     *
     * Now, the i-th step goes from Y[i] to Y[i + 1] */
    double curr_omega[3][3];
    double curr_ijk_ijk[3][3][3];
    double curr_ijk_xyz[3][3][3];
    double dts[3] = {0.0, dt / 2, dt};
    curr_ijk_xyz[0][0][0] = *ix;
    curr_ijk_xyz[0][0][1] = *iy;
    curr_ijk_xyz[0][0][2] = *iz;
    curr_ijk_xyz[0][1][0] = *jx;
    curr_ijk_xyz[0][1][1] = *jy;
    curr_ijk_xyz[0][1][2] = *jz;
    curr_ijk_xyz[0][2][0] = *kx;
    curr_ijk_xyz[0][2][1] = *ky;
    curr_ijk_xyz[0][2][2] = *kz;
    curr_ijk_ijk[0][0][0] = 1.0;
    curr_ijk_ijk[0][0][1] = 0.0;
    curr_ijk_ijk[0][0][2] = 0.0;
    curr_ijk_ijk[0][1][0] = 0.0;
    curr_ijk_ijk[0][1][1] = 1.0;
    curr_ijk_ijk[0][1][2] = 0.0;
    curr_ijk_ijk[0][2][0] = 0.0;
    curr_ijk_ijk[0][2][1] = 0.0;
    curr_ijk_ijk[0][2][2] = 1.0;
    curr_omega[0][0] = *omega**si;
    curr_omega[0][1] = *omega**sj;
    curr_omega[0][2] = *omega**sk;

    for (int curr = 0; curr < 2; curr++) {
        /******* CALCULATE all dY/dt (Y = each variable) ***************************/
        rebx_calc_torques(sim,index,M_ijk,I_ijk,curr_ijk_xyz[curr], dts[curr], dt);
        rebx_domega_dt(curr_omega[curr],M_ijk,I_ijk,domega_dts);
        rebx_dijk_dt(curr_ijk_ijk[curr],curr_omega[curr],dijk_dts);

        /******* CALCULATE all dYs ***************************/
        for (int i = 0; i < 3; i++){
            domega[i] = domega_dts[i] * dts[curr + 1];
            for (int j = 0; j < 3; j++){
                dijk[i][j] = dijk_dts[i][j] * dts[curr + 1];
            }
        }

        /******* CALCULATE all new Y's, store to curr_Y[curr + 1] ***********************/

        for (int i = 0; i < 3; i++) {
            curr_omega[curr + 1][i] = curr_omega[0][i] + domega[i];
            for (int j = 0; j < 3; j++) {
                curr_ijk_ijk[curr + 1][i][j] = curr_ijk_ijk[0][i][j] + dijk[i][j];
            }
            rebx_ijk_to_xyz(curr_ijk_ijk[curr + 1][i],curr_ijk_xyz[curr + 1][i],ijk_xyz);
        }
        renorm_mat(curr_ijk_xyz[curr + 1]);
    }

    /******* UPDATE ALL Y's ***************************/
    *omega = mag_vec(curr_omega[2]);
    renorm_vec(curr_omega[2]);
    *si = curr_omega[2][0];
    *sj = curr_omega[2][1];
    *sk = curr_omega[2][2];
    *ix = curr_ijk_xyz[2][0][0];
    *iy = curr_ijk_xyz[2][0][1];
    *iz = curr_ijk_xyz[2][0][2];
    *jx = curr_ijk_xyz[2][1][0];
    *jy = curr_ijk_xyz[2][1][1];
    *jz = curr_ijk_xyz[2][1][2];
    *kx = curr_ijk_xyz[2][2][0];
    *ky = curr_ijk_xyz[2][2][1];
    *kz = curr_ijk_xyz[2][2][2];
}

// runs checks on parameters. Returns 1 if error, 0 otherwise.
static int rebx_validate_params(struct reb_simulation* const sim, double* const ix, double* const iy, double* const iz, double* const jx, double* const jy,
    double* const jz, double* const kx, double* const ky, double* const kz, double* const si,
    double* const sj, double* const sk) {

    double tolerance = 1.e-15;

    double i_diff = fabs(*ix**ix + *iy**iy + *iz**iz - 1);
    double j_diff = fabs(*jx**jx + *jy**jy + *jz**jz - 1);
    double k_diff = fabs(*kx**kx + *ky**ky + *kz**kz - 1);
    double s_diff = fabs(*si**si + *sj**sj + *sk**sk - 1);

    double i_dot_j = *ix**jx + *iy**jy + *iz**jz;
    double i_dot_k = *ix**kx + *iy**ky + *iz**kz;
    double j_dot_k = *jx**kx + *jy**ky + *jz**kz;

    double i_cross_j_x = *iy**jz - *iz**jy;
    double i_cross_j_y = *iz**jx - *ix**jz;
    double i_cross_j_z = *ix**jy - *iy**jx;

    if (i_diff > tolerance || j_diff > tolerance || k_diff > tolerance || s_diff > tolerance) {
        reb_error(sim, "REBOUNDx Error: triaxial_torque: Vectors i, j, k, and s must be unit vectors. \n");
        return 1;
    }
    if (i_dot_j != 0 || i_dot_k != 0 || j_dot_k != 0) {
        reb_error(sim, "REBOUNDx Error: triaxial_torque: Vectors i, j, and k must be mutually orthogonal. \n");
        return 1;
    }
    if (i_cross_j_x != *kx || i_cross_j_y != *ky || i_cross_j_z != *kz){
        reb_error(sim, "REBOUNDx Error: triaxial_torque: The cross-product of vectors i and j must equal vector k. \n");
        return 1;
    }
}

void rebx_triaxial_torque(struct reb_simulation* const sim, struct rebx_operator* const triaxial_torque, const double dt){
    const int _N_real = sim->N - sim->N_var;
    for(int i=0; i<_N_real; i++){
        struct reb_particle* const p = &sim->particles[i];

        // check required params
        const double* const Ii = rebx_get_param(sim->extras, p->ap, "tt_Ii");
        if (Ii == NULL) {
            continue;
        }
        const double* const Ij = rebx_get_param(sim->extras, p->ap, "tt_Ij");
        if (Ij == NULL) {
            continue;
        }
        const double* const Ik = rebx_get_param(sim->extras, p->ap, "tt_Ik");
        if (Ik == NULL) {
            continue;
        }
        double* const omega = rebx_get_param(sim->extras, p->ap, "tt_omega");
        if (omega == NULL) {
            continue;
        }
        double* const ix = rebx_get_param(sim->extras, p->ap, "tt_ix");
        if (ix == NULL) {
            continue;
        }
        double* const iy = rebx_get_param(sim->extras, p->ap, "tt_iy");
        if (iy == NULL) {
            continue;
        }
        double* const iz = rebx_get_param(sim->extras, p->ap, "tt_iz");
        if (iz == NULL) {
            continue;
        }
        double* const jx = rebx_get_param(sim->extras, p->ap, "tt_jx");
        if (jx == NULL) {
            continue;
        }
        double* const jy = rebx_get_param(sim->extras, p->ap, "tt_jy");
        if (jy == NULL) {
            continue;
        }
        double* const jz = rebx_get_param(sim->extras, p->ap, "tt_jz");
        if (jz == NULL) {
            continue;
        }
        double* const kx = rebx_get_param(sim->extras, p->ap, "tt_kx");
        if (kx == NULL) {
            continue;
        }
        double* const ky = rebx_get_param(sim->extras, p->ap, "tt_ky");
        if (ky == NULL) {
            continue;
        }
        double* const kz = rebx_get_param(sim->extras, p->ap, "tt_kz");
        if (kz == NULL) {
            continue;
        }
        double* const si = rebx_get_param(sim->extras, p->ap, "tt_si");
        if (si == NULL) {
            continue;
        }
        double* const sj = rebx_get_param(sim->extras, p->ap, "tt_sj");
        if (sj == NULL) {
            continue;
        }
        double* const sk = rebx_get_param(sim->extras, p->ap, "tt_sk");
        if (sk == NULL) {
            continue;
        }

        // check validity of parameters if first timestep
        if (sim->t <= dt){
            if (rebx_validate_params(sim,ix,iy,iz,jx,jy,jz,kx,ky,kz,si,sj,sk) == 1) {
                return;
            }
        }

        /* rebx_update_spin_ijk(sim,1,i,ix,iy,iz,jx,jy,jz,kx,ky,kz,si,sj,sk,omega,*Ii,*Ij,*Ik,dt); */
        /* rebx_update_spin_ijk_euler(sim,i,ix,iy,iz,jx,jy,jz,kx,ky,kz,si,sj,sk,omega,*Ii,*Ij,*Ik,dt); */
        rebx_update_spin_ijk_mid(sim,i,ix,iy,iz,jx,jy,jz,kx,ky,kz,si,sj,sk,omega,*Ii,*Ij,*Ik,dt);
    }
}
