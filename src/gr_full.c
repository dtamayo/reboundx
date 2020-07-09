/**
 * @file    gr_full.c
 * @brief   Post-newtonian general relativity corrections for all bodies in the simulation
 * @author  Pengshuai (Sam) Shi, Hanno Rein, Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Pengshuai (Sam) Shi, Hanno Rein, Dan Tamayo
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
 * Authors                 P. Shi, H. Rein, D. Tamayo
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
 * Based on                `Newhall et al. 1983 <http://labs.adsabs.harvard.edu/adsabs/abs/1983A%26A...125..150N/>`_.
 * C Example               :ref:`c_example_gr`
 * Python Example          `GeneralRelativity.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GeneralRelativity.ipynb>`_.
 * ======================= ===============================================
 * 
 * This algorithm incorporates the first-order post-newtonian effects from all bodies in the system, and is necessary for multiple massive bodies like stellar binaries.
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
 * *None*
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_gr_full(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double C2, const double G, const int max_iterations, const int gravity_ignore_10){
    
    double a_const[N][3]; // array that stores the value of the constant term
    double a_newton[N][3]; // stores the Newtonian term
    double a_new[N][3]; // stores the newly calculated term
    double rs[N][N];
    double drs[N][N][3];

    for (int i=0; i<N; i++){
        // compute the Newtonian term 
        a_newton[i][0] = particles[i].ax;
        a_newton[i][1] = particles[i].ay;
        a_newton[i][2] = particles[i].az;
        a_new[i][0] = 0.;
        a_new[i][1] = 0.;
        a_new[i][2] = 0.;

        for(int j=i+1; j<N; j++){
            if (j!=i){
                drs[i][j][0] = particles[i].x - particles[j].x;
                drs[i][j][1] = particles[i].y - particles[j].y;
                drs[i][j][2] = particles[i].z - particles[j].z;
                drs[j][i][0] = -drs[i][j][0];
                drs[j][i][1] = -drs[i][j][1];
                drs[j][i][2] = -drs[i][j][2];
                rs[i][j] = sqrt(drs[i][j][0]*drs[i][j][0] + drs[i][j][1]*drs[i][j][1] + drs[i][j][2]*drs[i][j][2]);
                rs[j][i] = rs[i][j];
            }
        }
    }

    if (gravity_ignore_10){
        const double prefact = -G/(rs[0][1]*rs[0][1]*rs[0][1]);
        const double prefact0 = prefact*particles[0].m;
        const double prefact1 = prefact*particles[1].m;
        a_newton[0][0] += prefact1*drs[0][1][0];
        a_newton[0][1] += prefact1*drs[0][1][1];
        a_newton[0][2] += prefact1*drs[0][1][2];
        a_newton[1][0] -= prefact0*drs[0][1][0];
        a_newton[1][1] -= prefact0*drs[0][1][1];
        a_newton[1][2] -= prefact0*drs[0][1][2];
    }

    for (int i=0; i<N; i++){
        // then compute the constant terms:
        double a_constx = 0.;
        double a_consty = 0.;
        double a_constz = 0.;
        // 1st constant part
        for (int j = 0; j< N; j++){
            if (j != i){
                const double dxij = drs[i][j][0];
                const double dyij = drs[i][j][1];
                const double dzij = drs[i][j][2];
                const double rij2 = rs[i][j]*rs[i][j];
                const double rij3 = rij2*rs[i][j];
                
                double a1 = 0.;
                for (int k = 0; k< N; k++){
                    if (k != i){
                        a1 += (4./(C2)) * G*particles[k].m/rs[i][k];
                    }
                }

                double a2 = 0.;
                for (int l = 0; l< N; l++){
                    if (l != j){
                        a2 += (1./(C2)) * G*particles[l].m/rs[l][j];
                    }
                }

                double a3;
                double vi2 = particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy + particles[i].vz*particles[i].vz;
                a3 = -vi2/(C2);

                double a4;
                double vj2 = particles[j].vx*particles[j].vx + particles[j].vy*particles[j].vy + particles[j].vz*particles[j].vz;
                a4 = -2.*vj2/(C2);

                double a5;
                a5 = (4./(C2)) * (particles[i].vx*particles[j].vx + particles[i].vy*particles[j].vy + particles[i].vz*particles[j].vz); 
                
                double a6;
                double a6_0 = dxij*particles[j].vx + dyij*particles[j].vy + dzij*particles[j].vz;
                a6 = (3./(2.*C2)) * a6_0*a6_0/rij2;
                
                double factor1 = a1 + a2 + a3 + a4 + a5 + a6;
                 
                a_constx += G*particles[j].m*dxij*factor1/rij3;
                a_consty += G*particles[j].m*dyij*factor1/rij3;
                a_constz += G*particles[j].m*dzij*factor1/rij3;
        
                // 2nd constant part
                
                const double dvxij = particles[i].vx - particles[j].vx;
                const double dvyij = particles[i].vy - particles[j].vy;
                const double dvzij = particles[i].vz - particles[j].vz;
                    
                double factor2 = dxij*(4.*particles[i].vx-3.*particles[j].vx)+dyij*(4.*particles[i].vy-3.*particles[j].vy)+dzij*(4.*particles[i].vz-3.*particles[j].vz);

                a_constx += G*particles[j].m*factor2*dvxij/rij3/(C2);
                a_consty += G*particles[j].m*factor2*dvyij/rij3/(C2);
                a_constz += G*particles[j].m*factor2*dvzij/rij3/(C2);
            }
        }  

        a_const[i][0] = a_constx;
        a_const[i][1] = a_consty;
        a_const[i][2] = a_constz;
    }

    // Now running the substitution again and again through the loop below
    for (int k=0; k<10; k++){ // you can set k as how many substitution you want to make
        double a_old[N][3]; // initialize an arry that stores the information of previousu calculated accleration
        for (int i =0; i <N; i++){
            a_old[i][0] = a_new[i][0]; // when k = 0, a_new is the Newtownian term which calculated before
            a_old[i][1] = a_new[i][1];
            a_old[i][2] = a_new[i][2];
        }
        // now add on the non-constant term
        for (int i = 0; i < N; i++){ // a_j is used to update a_i and vice versa
            double non_constx = 0.;
            double non_consty = 0.;
            double non_constz = 0.;
            for (int j = 0; j < N; j++){
                if (j != i){
                    const double dxij = drs[i][j][0];
                    const double dyij = drs[i][j][1];
                    const double dzij = drs[i][j][2];
                    const double rij = rs[i][j];
                    const double rij3 = rij*rij*rij;
                    non_constx += (G*particles[j].m*dxij/rij3)*(dxij*(a_newton[j][0]+a_old[j][0])+dyij*(a_newton[j][1]+a_old[j][1])+\
                                dzij*(a_newton[j][2]+a_old[j][2]))/(2.*C2) + (7./(2.*C2))*G*particles[j].m*(a_newton[j][0]+a_old[j][0])/rij;
                    non_consty += (G*particles[j].m*dyij/rij3)*(dxij*(a_newton[j][0]+a_old[j][0])+dyij*(a_newton[j][1]+a_old[j][1])+\
                                dzij*(a_newton[j][2]+a_old[j][2]))/(2.*C2) + (7./(2.*C2))*G*particles[j].m*(a_newton[j][1]+a_old[j][1])/rij;
                    non_constz += (G*particles[j].m*dzij/rij3)*(dxij*(a_newton[j][0]+a_old[j][0])+dyij*(a_newton[j][1]+a_old[j][1])+\
                                dzij*(a_newton[j][2]+a_old[j][2]))/(2.*C2) + (7./(2.*C2))*G*particles[j].m*(a_newton[j][2]+a_old[j][2])/rij;
                }
            }
            a_new[i][0] = (a_const[i][0] + non_constx);
            a_new[i][1] = (a_const[i][1] + non_consty);
            a_new[i][2] = (a_const[i][2] + non_constz);
        }
        
        // break out loop if a_new is converging
        double maxdev = 0.;
        double dx, dy, dz;
        for (int i = 0; i < N; i++){
            dx = (fabs(a_new[i][0]) < 1.e-30) ? 0. : fabs(a_new[i][0] - a_old[i][0])/a_new[i][0];
            dy = (fabs(a_new[i][1]) < 1.e-30) ? 0. : fabs(a_new[i][1] - a_old[i][1])/a_new[i][1];
            dz = (fabs(a_new[i][2]) < 1.e-30) ? 0. : fabs(a_new[i][2] - a_old[i][2])/a_new[i][2];
            
            if (dx > maxdev) { maxdev = dx; }
            if (dy > maxdev) { maxdev = dy; }
            if (dz > maxdev) { maxdev = dz; }
        }
        
        if (maxdev < 1.e-30){
            break;
        }
        if (k==9){
            reb_warning(sim, "10 loops in rebx_gr_full did not converge.\n");
        }
    }
    // update acceleration in particles
    for (int i = 0; i <N;i++){
        particles[i].ax += a_new[i][0];
        particles[i].ay += a_new[i][1];
        particles[i].az += a_new[i][2];
    }
}

void rebx_gr_full(struct reb_simulation* const sim, struct rebx_force* const gr_full, struct reb_particle* const particles, const int N){
    double* c = rebx_get_param(sim->extras, gr_full->ap, "c");
    if (c == NULL){
        reb_error(sim, "REBOUNDx Error: Need to set speed of light in gr effect.  See examples in documentation.\n");
        return;
    }
    const double C2 = (*c)*(*c);
    const unsigned int gravity_ignore_10 = sim->gravity_ignore_terms==1;
    int* max_iterations = rebx_get_param(sim->extras, gr_full->ap, "max_iterations");
    if(max_iterations != NULL){
        rebx_calculate_gr_full(sim, particles, N, C2, sim->G, *max_iterations, gravity_ignore_10);
    }
    else{
        const int default_max_iterations = 10;
        rebx_calculate_gr_full(sim, particles, N, C2, sim->G, default_max_iterations, gravity_ignore_10);
    }
}

double rebx_gr_full_hamiltonian(struct rebx_extras* const rebx, const struct rebx_force* const force){
    if (rebx->sim == NULL){
        rebx_error(rebx, ""); // rebx_error gives meaningful err
        return 0;
    }
    struct reb_simulation* sim = rebx->sim;
    double* c = rebx_get_param(rebx, force->ap, "c");
    if (c == NULL){
        reb_error(sim, "Need to set speed of light in gr effect.  See examples in documentation.\n");
    }
    const double C2 = (*c)*(*c);
    const int N = sim->N - sim->N_var;
    const double G = sim->G;
    struct reb_particle* const particles = sim->particles;

	double e_kin = 0.;
	double e_pot = 0.;
	double e_pn  = 0.;
	
    struct reb_vec3d* vtilde = malloc(N*sizeof(*vtilde));
    struct reb_vec3d* vtilde_old = malloc(N*sizeof(*vtilde_old));
    for (int j=0; j<N;j++){
        vtilde[j].x = particles[j].vx;
        vtilde[j].y = particles[j].vy;
        vtilde[j].z = particles[j].vz;
    }
    
    for (int q=0; q<10;q++){
        for (int i=0;i<N;i++){
            struct reb_particle pi = particles[i];
            
            double vtildei2 = vtilde[i].x*vtilde[i].x + vtilde[i].y*vtilde[i].y + vtilde[i].z*vtilde[i].z;
            double A = (1. - 0.5*vtildei2/C2);
            
            double sumk = 0.;
            for (int k=0;k<N;k++){
                if (k!=i){
                    struct reb_particle pk = particles[k];
                    double xik = pk.x - pi.x;
                    double yik = pk.y - pi.y;
                    double zik = pk.z - pi.z;
                    double rik = sqrt(xik*xik + yik*yik + zik*zik);
                    sumk -= 2.*G*pk.m/rik;
                }
            }
            
            struct reb_vec3d dv_pn = {0.};
            for (int j=0;j<N;j++){
                if (j != i){
                    struct reb_particle pj = particles[j];
                    double xij = pj.x - pi.x;
                    double yij = pj.y - pi.y;
                    double zij = pj.z - pi.z;
                    double rij2 = xij*xij + yij*yij + zij*zij;
                    double rij = sqrt(rij2);
                    double rijdotvj = vtilde[j].x*xij + vtilde[j].y*yij + vtilde[j].z*zij;
                    double pfac = pj.m/rij;

                    dv_pn.x += pfac*(6.*vtilde[i].x - 7.*vtilde[j].x - rijdotvj*xij/rij2);
                    dv_pn.y += pfac*(6.*vtilde[i].y - 7.*vtilde[j].y - rijdotvj*yij/rij2);
                    dv_pn.z += pfac*(6.*vtilde[i].z - 7.*vtilde[j].z - rijdotvj*zij/rij2);
                }
            }

            dv_pn.x *= G/(2.*C2);
            dv_pn.y *= G/(2.*C2);
            dv_pn.z *= G/(2.*C2);

            vtilde[i].x = (pi.vx + dv_pn.x)/A;
            vtilde[i].y = (pi.vy + dv_pn.y)/A;
            vtilde[i].z = (pi.vz + dv_pn.z)/A;
        }
    }
    
    for (int i=0; i<N; i++){
        struct reb_particle p = particles[i];
        double vtildei2 = vtilde[i].x*vtilde[i].x + vtilde[i].y*vtilde[i].y + vtilde[i].z*vtilde[i].z;
        e_kin += 0.5*p.m*vtildei2;
    }

    for (int i=0;i<N;i++){
        struct reb_particle pi = particles[i];
        double sumk = 0.;
        for (int k=0;k<N;k++){
            if (k!=i){
                struct reb_particle pk = particles[k];
                double xik = pk.x - pi.x;
                double yik = pk.y - pi.y;
                double zik = pk.z - pi.z;
                double rik = sqrt(xik*xik + yik*yik + zik*zik);
                sumk -= 2.*G*pk.m/rik;
            }
        }
        
        double vtildei2 = vtilde[i].x*vtilde[i].x + vtilde[i].y*vtilde[i].y + vtilde[i].z*vtilde[i].z;

        for (int j=0;j<N;j++){
            if (j != i){
                struct reb_particle pj = particles[j];
                double xij = pj.x - pi.x;
                double yij = pj.y - pi.y;
                double zij = pj.z - pi.z;
                double rij2 = xij*xij + yij*yij + zij*zij;
                double rij = sqrt(rij2);
                double rijdotvj = vtilde[j].x*xij + vtilde[j].y*yij + vtilde[j].z*zij;
                double rijdotvi = vtilde[i].x*xij + vtilde[i].y*yij + vtilde[i].z*zij;
                double vidotvj = vtilde[i].x*vtilde[j].x + vtilde[i].y*vtilde[j].y + vtilde[i].z*vtilde[j].z;
                
                e_pn -= G/(4.*C2)*pi.m*pj.m/rij*(6.*vtildei2 - 7*vidotvj - rijdotvi*rijdotvj/rij2 + sumk);
            }
        }
        
        e_pn -= pi.m/(8.*C2)*vtildei2*vtildei2;
        
        for (int j=i+1; j<N; j++){ // classic full
            struct reb_particle pj = particles[j];
            double dx = pi.x - pj.x;
            double dy = pi.y - pj.y;
            double dz = pi.z - pj.z;	
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            
            e_pot -= G*pi.m*pj.m/r;
        }
    }
    
	return e_kin + e_pot + e_pn;
}
