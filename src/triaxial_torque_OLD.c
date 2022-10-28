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
 * Parameter reuqirements:
 * ============================
 * i, j, k, and s must be unit vectors (e.g. such that ix^2 + iy^2 + iz^2 = 1)
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

// computes time-derivative of spin vector omega_i,omega_j,omega_k using torque vector Mi, Mj, Mk according to Euler's equations
static void rebx_domega_dt(double omega_i, double omega_j, double omega_k, double* M_ijk, double Ii, double Ij,
    double Ik, double* domega_dts){

    domega_dts[0] = (M_ijk[0] + (Ij-Ik)*omega_j*omega_k) / Ii;
    domega_dts[1] = (M_ijk[1] + (Ik-Ii)*omega_k*omega_i) / Ij;
    domega_dts[2] = (M_ijk[2] + (Ii-Ij)*omega_i*omega_j) / Ik;
}

static void rebx_calc_torques_old(struct reb_simulation* const sim, int index, double* M_ijk, double Ii, double Ij, double Ik, double ix, 
    double iy, double iz, double jx, double jy, double jz, double kx, double ky, double kz){
    
    struct reb_particle* p = &sim->particles[index];
    struct reb_particle* torquer;
    double rx;
    double ry;
    double rz;
    double r;
    double r_dot_i;
    double r_dot_j;
    double r_dot_k;
    double prefac;

    const int _N_real = sim->N - sim->N_var;
	for(int i=0; i<_N_real; i++){
        if (i == index) {
            continue;
        }
        
        torquer = &sim->particles[i];
        rx = p->x - torquer->x;
        ry = p->y - torquer->y;
        rz = p->z - torquer->z;
        r = sqrt(rx*rx + ry*ry + rz*rz);
        prefac = 3 * sim->G * torquer->m / pow(r,5);

        r_dot_i = rx*ix + ry*iy + rz*iz;
        r_dot_j = rx*jx + ry*jy + rz*jz;
        r_dot_k = rx*kx + ry*ky + rz*kz;

        M_ijk[0] += prefac*(Ik-Ij)*r_dot_j*r_dot_k;
        M_ijk[1] += prefac*(Ii-Ik)*r_dot_k*r_dot_i;
        M_ijk[2] += prefac*(Ij-Ii)*r_dot_i*r_dot_j;
    }
}

/* updates spin direction vector and spin rate according to Euler's equations and torques applied by other bodies */
static void rebx_update_spin(struct reb_simulation* const sim, int index, double* const ix, double* const iy, double* const iz, 
    double* const jx, double* const jy, double* const jz, double* const kx, double* const ky, double* const kz, 
    double* const si, double* const sj, double* const sk, double* const omega, const double Ii, const double Ij, 
    const double Ik, const double dt) {
    
    // calculate torques in i,j,k 
    double M_ijk[3];
    M_ijk[0] = 0.0;
    M_ijk[1] = 0.0;
    M_ijk[2] = 0.0;
    // rebx_calc_torques_old(sim,index,M_ijk,Ii,Ij,Ik,*ix,*iy,*iz,*jx,*jy,*jz,*kx,*ky,*kz);

    // change sx sy sz to ijk basis
    // double det_ijk = *ix*(*jy**kz-*ky**jz) + *jx*(*ky**iz-*iy**kz) + *kx*(*iy**jz-*jy**iz); // determinant of ijk unit vector matrix
    
    // double inv_ijk[3][3]; // inverse of ijk unit vector matrix
    // inv_ijk[0][0] = (*jy**kz-*jz**ky) / det_ijk;
    // inv_ijk[0][1] = (*jz**kx-*jx**kz) / det_ijk;
    // inv_ijk[0][2] = (*jx**ky-*jy**kx) / det_ijk;
    // inv_ijk[1][0] = (*iz**ky-*iy**kz) / det_ijk;
    // inv_ijk[1][1] = (*ix**kz-*iz**kx) / det_ijk;
    // inv_ijk[1][2] = (*iy**kx-*ix**ky) / det_ijk;
    // inv_ijk[2][0] = (*iy**jz-*jy**iz) / det_ijk;
    // inv_ijk[2][1] = (*jx**iz-*jz**ix) / det_ijk;
    // inv_ijk[2][2] = (*jy**ix-*iy**jx) / det_ijk;

    // double omega_i = *omega*(*sx*inv_ijk[0][0] + *sy*inv_ijk[0][1] + *sz*inv_ijk[0][2]);
    // double omega_j = *omega*(*sx*inv_ijk[1][0] + *sy*inv_ijk[1][1] + *sz*inv_ijk[1][2]);
    // double omega_k = *omega*(*sx*inv_ijk[2][0] + *sy*inv_ijk[2][1] + *sz*inv_ijk[2][2]);

    double omega_i = *omega**si;
    double omega_j = *omega**sj;
    double omega_k = *omega**sk;

    /* matrix for all calculations of slope
    first dimension: which RK derivation calculation (1-4) 
    dsecond dimension: component of omega vector (i,j,k) */
    double domega_dts[4][3];

    // Four Runge Kutta calculations
    rebx_domega_dt(omega_i,omega_j,omega_k,M_ijk,Ii,Ij,Ik,domega_dts[0]);
    rebx_domega_dt(omega_i + (domega_dts[0][0] * dt * 0.5),omega_j + (domega_dts[0][1] * dt * 0.5),omega_k + (domega_dts[0][2] * dt * 0.5),M_ijk,Ii,Ij,Ik,domega_dts[1]);
    rebx_domega_dt(omega_i + (domega_dts[1][0] * dt * 0.5),omega_j + (domega_dts[1][1] * dt * 0.5),omega_k + (domega_dts[1][2] * dt * 0.5),M_ijk,Ii,Ij,Ik,domega_dts[2]);
    rebx_domega_dt(omega_i + (domega_dts[2][0] * dt),omega_j + (domega_dts[2][1] * dt),omega_k + (domega_dts[2][2] * dt),M_ijk,Ii,Ij,Ik,domega_dts[3]);

    // calculate domega
    double domega[3];
    for (int i = 0; i < 3; i++){
        domega[i] = (domega_dts[0][i] + 2*domega_dts[1][i] + 2*domega_dts[2][i] + domega_dts[3][i]) * dt / 6;
    }

    omega_i += domega[0];
    omega_j += domega[1];
    omega_k += domega[2];

    *omega = sqrt(omega_i*omega_i + omega_j*omega_j + omega_k*omega_k);
    *si = omega_i / *omega;
    *sj = omega_j / *omega;
    *sk = omega_k / *omega;

    // double si = omega_i / *omega;
    // double sj = omega_j / *omega;
    // double sk = omega_k / *omega;

    // convert back to xyz, update sx, sy, sz
    // *sx = *ix*si + *jx*sj + *kx*sk;
    // *sy = *iy*si + *jy*sj + *ky*sk;
    // *sz = *iz*si + *jz*sj + *kz*sk;
}

// computes time-derivative of vector i,j,k
static void rebx_dijk_dt(double i, double j, double k, double si, double sj, double sk, double omega, double* dijk_dt){
    dijk_dt[0] = omega*(sj*k - sk*j);
    dijk_dt[1] = omega*(sk*i - si*k);
    dijk_dt[2] = omega*(si*j - sj*i);
}

// updates i, j, k vectors based on spin rate omega, spin vector s, and timestep dt
static void rebx_update_ijk(double* const ix, double* const iy, double* const iz, double* const jx, double* const jy, 
    double* const jz, double* const kx, double* const ky, double* const kz, double* const si, double* const sj, 
    double* const sk, double* const omega, const double dt) {

    /* matrix for all calculations of slope
    first dimension: which RK derivation calculation (1-4)
    second dimension: vector (i_hat,j_hat,k_hat)
    third dimension: component of vector (i,j,k) */
    double dijk_dts[4][3][3];

    // first Runge Kutta calculation
    rebx_dijk_dt(1.0,0.0,0.0,*si,*sj,*sk,*omega,dijk_dts[0][0]); // i_hat
    rebx_dijk_dt(0.0,1.0,0.0,*si,*sj,*sk,*omega,dijk_dts[0][1]); // j_hat
    rebx_dijk_dt(0.0,0.0,1.0,*si,*sj,*sk,*omega,dijk_dts[0][2]); // k_hat

    // second Runge Kutta calculation
    rebx_dijk_dt(1.0+(dijk_dts[0][0][0]*dt*0.5),0.0+(dijk_dts[0][0][1]*dt*0.5),0.0+(dijk_dts[0][0][2]*dt*0.5),*si,*sj,*sk,*omega,dijk_dts[1][0]); // i_hat
    rebx_dijk_dt(0.0+(dijk_dts[0][1][0]*dt*0.5),1.0+(dijk_dts[0][1][1]*dt*0.5),0.0+(dijk_dts[0][1][2]*dt*0.5),*si,*sj,*sk,*omega,dijk_dts[1][1]); // j_hat
    rebx_dijk_dt(0.0+(dijk_dts[0][2][0]*dt*0.5),0.0+(dijk_dts[0][2][1]*dt*0.5),1.0+(dijk_dts[0][2][2]*dt*0.5),*si,*sj,*sk,*omega,dijk_dts[1][2]); // k_hat

    // third Runge Kutta calculation
    rebx_dijk_dt(1.0+(dijk_dts[1][0][0]*dt*0.5),0.0+(dijk_dts[1][0][1]*dt*0.5),0.0+(dijk_dts[1][0][2]*dt*0.5),*si,*sj,*sk,*omega,dijk_dts[2][0]); // i_hat
    rebx_dijk_dt(0.0+(dijk_dts[1][1][0]*dt*0.5),1.0+(dijk_dts[1][1][1]*dt*0.5),0.0+(dijk_dts[1][1][2]*dt*0.5),*si,*sj,*sk,*omega,dijk_dts[2][1]); // j_hat
    rebx_dijk_dt(0.0+(dijk_dts[1][2][0]*dt*0.5),0.0+(dijk_dts[1][2][1]*dt*0.5),1.0+(dijk_dts[1][2][2]*dt*0.5),*si,*sj,*sk,*omega,dijk_dts[2][2]); // k_hat
    
    // fourth Runge Kutta calculation
    rebx_dijk_dt(1.0+(dijk_dts[2][0][0]*dt),0.0+(dijk_dts[2][0][1]*dt),0.0+(dijk_dts[2][0][2]*dt),*si,*sj,*sk,*omega,dijk_dts[3][0]); // i_hat
    rebx_dijk_dt(0.0+(dijk_dts[2][1][0]*dt),1.0+(dijk_dts[2][1][1]*dt),0.0+(dijk_dts[2][1][2]*dt),*si,*sj,*sk,*omega,dijk_dts[3][1]); // j_hat
    rebx_dijk_dt(0.0+(dijk_dts[2][2][0]*dt),0.0+(dijk_dts[2][2][1]*dt),1.0+(dijk_dts[2][2][2]*dt),*si,*sj,*sk,*omega,dijk_dts[3][2]); // k_hat

    // calculate dijk
    double dijks[3][3];
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            dijks[i][j] = (dijk_dts[0][i][j] + 2*dijk_dts[1][i][j] + 2*dijk_dts[2][i][j] + dijk_dts[3][i][j]) * dt / 6;
        }
    }

    double ix_ = *ix + dijks[0][0]**ix + dijks[0][1]**jx + dijks[0][2]**kx;
    double iy_ = *iy + dijks[0][0]**iy + dijks[0][1]**jy + dijks[0][2]**ky;
    double iz_ = *iz + dijks[0][0]**iz + dijks[0][1]**jz + dijks[0][2]**kz;
    double jx_ = *jx + dijks[1][0]**ix + dijks[1][1]**jx + dijks[1][2]**kx;
    double jy_ = *jy + dijks[1][0]**iy + dijks[1][1]**jy + dijks[1][2]**ky;
    double jz_ = *jz + dijks[1][0]**iz + dijks[1][1]**jz + dijks[1][2]**kz;
    double kx_ = *kx + dijks[2][0]**ix + dijks[2][1]**jx + dijks[2][2]**kx;
    double ky_ = *ky + dijks[2][0]**iy + dijks[2][1]**jy + dijks[2][2]**ky;
    double kz_ = *kz + dijks[2][0]**iz + dijks[2][1]**jz + dijks[2][2]**kz;

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

    /*********************************************************/

    // double dtheta = *omega * dt;
    // double sin_theta = sin(dtheta);
    // double cos_theta = sqrt(1-(sin_theta*sin_theta));

    // double r[3][3]; // rotation matrix ([row][col])
    // r[0][0] = cos_theta + (*si * *si * (1-cos_theta));
    // r[0][1] = (*si * *sj * (1-cos_theta)) - (*sk * sin_theta);
    // r[0][2] = (*si * *sk * (1-cos_theta)) + (*sj * sin_theta);
    // r[1][0] = (*si * *sj * (1-cos_theta)) + (*sk * sin_theta);
    // r[1][1] = cos_theta + (*sj * *sj * (1-cos_theta));
    // r[1][2] = (*sk * *sj * (1-cos_theta)) - (*si * sin_theta);
    // r[2][0] = (*si * *sk * (1-cos_theta)) - (*sj * sin_theta);
    // r[2][1] = (*sk * *sj * (1-cos_theta)) + (*si * sin_theta);
    // r[2][2] = cos_theta + (*sk * *sk * (1-cos_theta));

    // double ii = r[0][0];
    // double ij = r[1][0];
    // double ik = r[2][0];
    // double ji = r[0][1];
    // double jj = r[1][1];
    // double jk = r[2][1];
    // double ki = r[0][2];
    // double kj = r[1][2];
    // double kk = r[2][2];

    // double ix_ = ii**ix + ij**jx + ik**kx;
    // double iy_ = ii**iy + ij**jy + ik**ky;
    // double iz_ = ii**iz + ij**jz + ik**kz;
    // double jx_ = ji**ix + jj**jx + jk**kx;
    // double jy_ = ji**iy + jj**jy + jk**ky;
    // double jz_ = ji**iz + jj**jz + jk**kz;
    // double kx_ = ki**ix + kj**jx + kk**kx;
    // double ky_ = ki**iy + kj**jy + kk**ky;
    // double kz_ = ki**iz + kj**jz + kk**kz;

    // *ix = ix_;
    // *iy = iy_;
    // *iz = iz_;
    // *jx = jx_;
    // *jy = jy_;
    // *jz = jz_;
    // *kx = kx_;
    // *ky = ky_;
    // *kz = kz_;

    /*********************************************************/

    // convert to sx, sy, sz
    // double sx = *si**ix + *sj**jx + *sk**kx;
    // double sy = *si**iy + *sj**jy + *sk**ky;
    // double sz = *si**iz + *sj**jz + *sk**kz;

    // double r[3][3]; // rotation matrix ([row][col])
    // r[0][0] = cos_theta + (*sx * *sx * (1-cos_theta));
    // r[0][1] = (*sx * *sy * (1-cos_theta)) - (*sz * sin_theta);
    // r[0][2] = (*sx * *sz * (1-cos_theta)) + (*sy * sin_theta);
    // r[1][0] = (*sx * *sy * (1-cos_theta)) + (*sz * sin_theta);
    // r[1][1] = cos_theta + (*sy * *sy * (1-cos_theta));
    // r[1][2] = (*sz * *sy * (1-cos_theta)) - (*sx * sin_theta);
    // r[2][0] = (*sx * *sz * (1-cos_theta)) - (*sy * sin_theta);
    // r[2][1] = (*sz * *sy * (1-cos_theta)) + (*sx * sin_theta);
    // r[2][2] = cos_theta + (*sz * *sz * (1-cos_theta));

    // double ix_ = *ix * r[0][0] + *iy * r[0][1] + *iz * r[0][2];
    // double iy_ = *ix * r[1][0] + *iy * r[1][1] + *iz * r[1][2];
    // double iz_ = *ix * r[2][0] + *iy * r[2][1] + *iz * r[2][2];
    // double jx_ = *jx * r[0][0] + *jy * r[0][1] + *jz * r[0][2];
    // double jy_ = *jx * r[1][0] + *jy * r[1][1] + *jz * r[1][2];
    // double jz_ = *jx * r[2][0] + *jy * r[2][1] + *jz * r[2][2];
    // double kx_ = *kx * r[0][0] + *ky * r[0][1] + *kz * r[0][2];
    // double ky_ = *kx * r[1][0] + *ky * r[1][1] + *kz * r[1][2];
    // double kz_ = *kx * r[2][0] + *ky * r[2][1] + *kz * r[2][2];

    // *ix = ix_;
    // *iy = iy_;
    // *iz = iz_;
    // *jx = jx_;
    // *jy = jy_;
    // *jz = jz_;
    // *kx = kx_;
    // *ky = ky_;
    // *kz = kz_;
}

static void rebx_calc_torques(struct reb_simulation* const sim, int index, double* xyz, double* M_ijk, double* I_ijk,
    double* i_xyz, double* j_xyz, double* k_xyz){
    
    struct reb_particle* torquer;
    double rx;
    double ry;
    double rz;
    double r;
    double r_dot_i;
    double r_dot_j;
    double r_dot_k;
    double prefac;

    const int _N_real = sim->N - sim->N_var;
	for(int i=0; i<_N_real; i++){
        if (i == index) {
            continue;
        }
        
        torquer = &sim->particles[i];
        rx = xyz[0] - torquer->x;
        ry = xyz[1] - torquer->y;
        rz = xyz[2] - torquer->z;
        r = sqrt(rx*rx + ry*ry + rz*rz);
        prefac = 3 * sim->G * torquer->m / pow(r,5);

        r_dot_i = rx*i_xyz[0] + ry*i_xyz[1] + rz*i_xyz[2];
        r_dot_j = rx*j_xyz[0] + ry*j_xyz[1] + rz*j_xyz[2];
        r_dot_k = rx*k_xyz[0] + ry*k_xyz[1] + rz*k_xyz[2];

        M_ijk[0] += prefac*(I_ijk[2]-I_ijk[1])*r_dot_j*r_dot_k;
        M_ijk[1] += prefac*(I_ijk[0]-I_ijk[2])*r_dot_k*r_dot_i;
        M_ijk[2] += prefac*(I_ijk[1]-I_ijk[0])*r_dot_i*r_dot_j;
    }
}

// updates spin vector, omega, and ijk in lockstep using 4th order Runge Kutta
// static void rebx_update_spin_ijk(struct reb_simulation* const sim, int index, double* const ix, double* const iy, double* const iz, 
//     double* const jx, double* const jy, double* const jz, double* const kx, double* const ky, double* const kz, double* const si, 
//     double* const sj, double* const sk, double* const omega, const double Ii, const double Ij, const double Ik, const double dt){

//     double I_ijk[3];
//     I_ijk[0] = Ii;
//     I_ijk[1] = Ij;
//     I_ijk[2] = Ik;

//     // initialize matrices for all R-K calculations
//     double rk_xyz[4][3]; // xyz of particle
//     double rk_M_ijk[4][3]; // ijk components of torque on body
//     double rk_omega_ijk[4][3]; // ijk components of spin vector (total omega vector)
//     double rk_ijk_xyz[4][3][3]; // xyz components of each ijk vector

//     struct reb_particle* p = &sim->particles[index];
//     rk_xyz[0][0] = p->x;
//     rk_xyz[0][1] = p->y;
//     rk_xyz[0][2] = p->z;

//     rk_omega_ijk[0][0] = *omega**si;
//     rk_omega_ijk[0][1] = *omega**sj;
//     rk_omega_ijk[0][2] = *omega**sk;

//     rk_ijk_xyz[0][0][0] = *ix;
//     rk_ijk_xyz[0][0][1] = *iy;
//     rk_ijk_xyz[0][0][2] = *iz;

//     rk_ijk_xyz[0][1][0] = *jx;
//     rk_ijk_xyz[0][1][1] = *jy;
//     rk_ijk_xyz[0][1][2] = *jz;

//     rk_ijk_xyz[0][2][0] = *kx;
//     rk_ijk_xyz[0][2][1] = *ky;
//     rk_ijk_xyz[0][2][2] = *kz;

//     for (int i=0; i < 4; i++) {
//         for (int j=0; j < 3; j++) {
//             M_ijk[i][j] = 0.0;
//         }
//     }

//     // First lockstep Runge-Kutta calculations
//     rebx_calc_torques(sim,index,rk_xyz[0],M_ijk[0],I_ijk,rk_ijk_xyz[0][0],rk_ijk_xyz[0][1],rk_ijk_xyz[0][2]);

//     // Continue editing here

//     /* matrix for all calculations of slope
//     first dimension: which RK derivation calculation (1-4) 
//     dsecond dimension: component of omega vector (i,j,k) */
//     double domega_dts[4][3];

//     // Four Runge Kutta calculations
//     rebx_domega_dt(omega_i,omega_j,omega_k,M_ijk,Ii,Ij,Ik,domega_dts[0]);
//     rebx_domega_dt(omega_i + (domega_dts[0][0] * dt * 0.5),omega_j + (domega_dts[0][1] * dt * 0.5),omega_k + (domega_dts[0][2] * dt * 0.5),M_ijk,Ii,Ij,Ik,domega_dts[1]);
//     rebx_domega_dt(omega_i + (domega_dts[1][0] * dt * 0.5),omega_j + (domega_dts[1][1] * dt * 0.5),omega_k + (domega_dts[1][2] * dt * 0.5),M_ijk,Ii,Ij,Ik,domega_dts[2]);
//     rebx_domega_dt(omega_i + (domega_dts[2][0] * dt),omega_j + (domega_dts[2][1] * dt),omega_k + (domega_dts[2][2] * dt),M_ijk,Ii,Ij,Ik,domega_dts[3]);

//     // calculate domega
//     double domega[3];
//     for (int i = 0; i < 3; i++){
//         domega[i] = (domega_dts[0][i] + 2*domega_dts[1][i] + 2*domega_dts[2][i] + domega_dts[3][i]) * dt / 6;
//     }

//     omega_i += domega[0];
//     omega_j += domega[1];
//     omega_k += domega[2];

//     *omega = sqrt(omega_i*omega_i + omega_j*omega_j + omega_k*omega_k);
//     *si = omega_i / *omega;
//     *sj = omega_j / *omega;
//     *sk = omega_k / *omega;

//     // double si = omega_i / *omega;
//     // double sj = omega_j / *omega;
//     // double sk = omega_k / *omega;

//     // convert back to xyz, update sx, sy, sz
//     // *sx = *ix*si + *jx*sj + *kx*sk;
//     // *sy = *iy*si + *jy*sj + *ky*sk;
//     // *sz = *iz*si + *jz*sj + *kz*sk;

//     /* matrix for all calculations of slope
//     first dimension: which RK derivation calculation (1-4)
//     second dimension: vector (i_hat,j_hat,k_hat)
//     third dimension: component of vector (i,j,k) */
//     double dijk_dts[4][3][3];

//     // first Runge Kutta calculation
//     rebx_dijk_dt(1.0,0.0,0.0,*si,*sj,*sk,*omega,dijk_dts[0][0]); // i_hat
//     rebx_dijk_dt(0.0,1.0,0.0,*si,*sj,*sk,*omega,dijk_dts[0][1]); // j_hat
//     rebx_dijk_dt(0.0,0.0,1.0,*si,*sj,*sk,*omega,dijk_dts[0][2]); // k_hat

//     // second Runge Kutta calculation
//     rebx_dijk_dt(1.0+(dijk_dts[0][0][0]*dt*0.5),0.0+(dijk_dts[0][0][1]*dt*0.5),0.0+(dijk_dts[0][0][2]*dt*0.5),*si,*sj,*sk,*omega,dijk_dts[1][0]); // i_hat
//     rebx_dijk_dt(0.0+(dijk_dts[0][1][0]*dt*0.5),1.0+(dijk_dts[0][1][1]*dt*0.5),0.0+(dijk_dts[0][1][2]*dt*0.5),*si,*sj,*sk,*omega,dijk_dts[1][1]); // j_hat
//     rebx_dijk_dt(0.0+(dijk_dts[0][2][0]*dt*0.5),0.0+(dijk_dts[0][2][1]*dt*0.5),1.0+(dijk_dts[0][2][2]*dt*0.5),*si,*sj,*sk,*omega,dijk_dts[1][2]); // k_hat

//     // third Runge Kutta calculation
//     rebx_dijk_dt(1.0+(dijk_dts[1][0][0]*dt*0.5),0.0+(dijk_dts[1][0][1]*dt*0.5),0.0+(dijk_dts[1][0][2]*dt*0.5),*si,*sj,*sk,*omega,dijk_dts[2][0]); // i_hat
//     rebx_dijk_dt(0.0+(dijk_dts[1][1][0]*dt*0.5),1.0+(dijk_dts[1][1][1]*dt*0.5),0.0+(dijk_dts[1][1][2]*dt*0.5),*si,*sj,*sk,*omega,dijk_dts[2][1]); // j_hat
//     rebx_dijk_dt(0.0+(dijk_dts[1][2][0]*dt*0.5),0.0+(dijk_dts[1][2][1]*dt*0.5),1.0+(dijk_dts[1][2][2]*dt*0.5),*si,*sj,*sk,*omega,dijk_dts[2][2]); // k_hat
    
//     // fourth Runge Kutta calculation
//     rebx_dijk_dt(1.0+(dijk_dts[2][0][0]*dt),0.0+(dijk_dts[2][0][1]*dt),0.0+(dijk_dts[2][0][2]*dt),*si,*sj,*sk,*omega,dijk_dts[3][0]); // i_hat
//     rebx_dijk_dt(0.0+(dijk_dts[2][1][0]*dt),1.0+(dijk_dts[2][1][1]*dt),0.0+(dijk_dts[2][1][2]*dt),*si,*sj,*sk,*omega,dijk_dts[3][1]); // j_hat
//     rebx_dijk_dt(0.0+(dijk_dts[2][2][0]*dt),0.0+(dijk_dts[2][2][1]*dt),1.0+(dijk_dts[2][2][2]*dt),*si,*sj,*sk,*omega,dijk_dts[3][2]); // k_hat

//     // calculate dijk
//     double dijks[3][3];
//     for (int i = 0; i < 3; i++){
//         for (int j = 0; j < 3; j++){
//             dijks[i][j] = (dijk_dts[0][i][j] + 2*dijk_dts[1][i][j] + 2*dijk_dts[2][i][j] + dijk_dts[3][i][j]) * dt / 6;
//         }
//     }

//     double ix_ = *ix + dijks[0][0]**ix + dijks[0][1]**jx + dijks[0][2]**kx;
//     double iy_ = *iy + dijks[0][0]**iy + dijks[0][1]**jy + dijks[0][2]**ky;
//     double iz_ = *iz + dijks[0][0]**iz + dijks[0][1]**jz + dijks[0][2]**kz;
//     double jx_ = *jx + dijks[1][0]**ix + dijks[1][1]**jx + dijks[1][2]**kx;
//     double jy_ = *jy + dijks[1][0]**iy + dijks[1][1]**jy + dijks[1][2]**ky;
//     double jz_ = *jz + dijks[1][0]**iz + dijks[1][1]**jz + dijks[1][2]**kz;
//     double kx_ = *kx + dijks[2][0]**ix + dijks[2][1]**jx + dijks[2][2]**kx;
//     double ky_ = *ky + dijks[2][0]**iy + dijks[2][1]**jy + dijks[2][2]**ky;
//     double kz_ = *kz + dijks[2][0]**iz + dijks[2][1]**jz + dijks[2][2]**kz;

//     // re-normalize
//     double i_mag = sqrt(ix_*ix_ + iy_*iy_ + iz_*iz_);
//     double j_mag = sqrt(jx_*jx_ + jy_*jy_ + jz_*jz_);
//     double k_mag = sqrt(kx_*kx_ + ky_*ky_ + kz_*kz_);

//     *ix = ix_ / i_mag;
//     *iy = iy_ / i_mag;
//     *iz = iz_ / i_mag;
//     *jx = jx_ / j_mag;
//     *jy = jy_ / j_mag;
//     *jz = jz_ / j_mag;
//     *kx = kx_ / k_mag;
//     *ky = ky_ / k_mag;
//     *kz = kz_ / k_mag;
// }

// runs checks on parameters. Returns 1 if error, 0 otherwise.
static int rebx_validate_params(struct reb_simulation* const sim, double* const ix, double* const iy, double* const iz, double* const jx, double* const jy,
    double* const jz, double* const kx, double* const ky, double* const kz, double* const si,
    double* const sj, double* const sk) {

    double tolerance = 1.e-15;

    double i_diff = abs(*ix**ix + *iy**iy + *iz**iz - 1);
    double j_diff = abs(*jx**jx + *jy**jy + *jz**jz - 1);
    double k_diff = abs(*kx**kx + *ky**ky + *kz**kz - 1);
    double s_diff = abs(*si**si + *sj**sj + *sk**sk - 1);

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

void rebx_triaxial_torque_OLD(struct reb_simulation* const sim, struct rebx_operator* const triaxial_torque, const double dt){
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
        
        // rebx_update_spin_ijk(sim,i,ix,iy,iz,jx,jy,jz,kx,ky,kz,si,sj,sk,omega,*Ii,*Ij,*Ik,dt);

        rebx_update_ijk(ix,iy,iz,jx,jy,jz,kx,ky,kz,si,sj,sk,omega,dt);
        rebx_update_spin(sim,i,ix,iy,iz,jx,jy,jz,kx,ky,kz,si,sj,sk,omega,*Ii,*Ij,*Ik,dt);
    }
}
