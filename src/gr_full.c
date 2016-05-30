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
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gr_full.h"
#include "rebound.h"
#include "reboundx.h"

struct rebx_params_gr_full* rebx_add_gr_full(struct rebx_extras* rebx, double c){
	struct rebx_params_gr_full* params = malloc(sizeof(*params));
	params->c = c;
    int force_is_velocity_dependent = 1;
    rebx_add_force(rebx, params, "gr_full", rebx_gr_full, force_is_velocity_dependent);
    return params;
}

void rebx_gr_full(struct reb_simulation* const sim, struct rebx_effect* const gr){
    const struct rebx_params_gr_full* const params = gr->paramsPtr;
    const double C = params->c;
    const int _N_real = sim->N - sim->N_var;
    const double G = sim->G;
    struct reb_particle* const particles = sim->particles;
    const unsigned int _gravity_ignore_10 = sim->gravity_ignore_10;

    double a_const[_N_real][3]; // array that stores the value of the constant term
    double a_newton[_N_real][3]; // stores the Newtonian term
    double a_new[_N_real][3]; // stores the newly calculated term

    for (int i=0; i<_N_real; i++){
        // compute the Newtonian term 
        a_newton[i][0] = particles[i].ax;
        a_newton[i][1] = particles[i].ay;
        a_newton[i][2] = particles[i].az;
        a_new[i][0] = a_newton[i][0]; // we want to use Newtonian term as our first substitution, hence the assignment here
        a_new[i][1] = a_newton[i][1];
        a_new[i][2] = a_newton[i][1];
    }
    if (_gravity_ignore_10){
        const double dx = particles[0].x - particles[1].x;
        const double dy = particles[0].y - particles[1].y;
        const double dz = particles[0].z - particles[1].z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
        const double prefact = -G/(r2*r);
        const double prefact0 = prefact*particles[0].m;
        const double prefact1 = prefact*particles[1].m;
        a_newton[0][0] += prefact1*dx;
        a_newton[0][1] += prefact1*dy;
        a_newton[0][2] += prefact1*dz;
        a_newton[1][0] -= prefact0*dx;
        a_newton[1][1] -= prefact0*dy;
        a_newton[1][2] -= prefact0*dz;

    }

    for (int i=0; i<_N_real; i++){
        // then compute the constant terms:
        double a_constx = 0.;
        double a_consty = 0.;
        double a_constz = 0.;
        // 1st constant part
        for (int j = 0; j< _N_real; j++){
            if (j != i){
                const double dxij = particles[i].x - particles[j].x;
                const double dyij = particles[i].y - particles[j].y;
                const double dzij = particles[i].z - particles[j].z;
                const double r2ij = dxij*dxij + dyij*dyij + dzij*dzij;
                const double rij = sqrt(r2ij);
                
                double a1 = 0.;
                for (int k = 0; k< _N_real; k++){
                    if (k != i){
                        const double dxik = particles[i].x - particles[k].x;
                        const double dyik = particles[i].y - particles[k].y;
                        const double dzik = particles[i].z - particles[k].z;
                        const double r2ik = dxik*dxik + dyik*dyik + dzik*dzik;
                        const double rik = sqrt(r2ik);
                        a1 += (4./(C*C)) * G*particles[k].m/rik;
                    }
                }

                double a2 = 0.;
                for (int l = 0; l< _N_real; l++){
                    if (l != j){
                        const double dxlj = particles[l].x - particles[j].x;
                        const double dylj = particles[l].y - particles[j].y;
                        const double dzlj = particles[l].z - particles[j].z;
                        const double r2lj = dxlj*dxlj + dylj*dylj + dzlj*dzlj;
                        const double rlj = sqrt(r2lj);
                        a2 += (1./(C*C)) * G*particles[l].m/rlj;
                    }
                }

                double a3;
                double vi2 = particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy + particles[i].vz*particles[i].vz;
                a3 = -vi2/(C*C);

                double a4;
                double vj2 = particles[j].vx*particles[j].vx + particles[j].vy*particles[j].vy + particles[j].vz*particles[j].vz;
                a4 = -2.*vj2/(C*C);

                double a5;
                a5 = (4./(C*C)) * (particles[i].vx*particles[j].vx + particles[i].vy*particles[j].vy + particles[i].vz*particles[j].vz); 
                
                double a6;
                double a6_0 = dxij*particles[j].vx + dyij*particles[j].vy + dzij*particles[j].vz;
                a6 = (3./(2.*C*C)) * a6_0*a6_0/r2ij;
                
                double factor1 = -1. + a1 + a2 + a3 + a4 + a5 + a6;
                 
                a_constx += G*particles[j].m*dxij*factor1/(r2ij*rij);
                a_consty += G*particles[j].m*dyij*factor1/(r2ij*rij);
                a_constz += G*particles[j].m*dzij*factor1/(r2ij*rij);
            }
        }   
        // 2nd consant part
        for (int j = 0; j< _N_real; j++){
            if (j != i){
                const double dxij = particles[i].x - particles[j].x;
                const double dyij = particles[i].y - particles[j].y;
                const double dzij = particles[i].z - particles[j].z;
                const double r2ij = dxij*dxij + dyij*dyij + dzij*dzij;
                const double rij = sqrt(r2ij);
                const double dvxij = particles[i].vx - particles[j].vx;
                const double dvyij = particles[i].vy - particles[j].vy;
                const double dvzij = particles[i].vz - particles[j].vz;
                    
                double factor2 = dxij*(4.*particles[i].vx-3.*particles[j].vx)+dyij*(4.*particles[i].vy-3.*particles[j].vy)+dzij*(4.*particles[i].vz-3.*particles[j].vz);

                a_constx += G*particles[j].m*factor2*dvxij/(r2ij*rij)/(C*C);
                a_consty += G*particles[j].m*factor2*dvyij/(r2ij*rij)/(C*C);
                a_constz += G*particles[j].m*factor2*dvzij/(r2ij*rij)/(C*C);
            }
        }
        a_const[i][0] = a_constx;
        a_const[i][1] = a_consty;
        a_const[i][2] = a_constz;
    }
    //fprintf(stderr,"%.16e\t%.16e\n",particles[1].ax, a_newton[1][0]);

    // Now running the substitution again and again through the loop below
    for (int k=0; k<10; k++){ // you can set k as how many substitution you want to make
        double a_old[_N_real][3]; // initialize an arry that stores the information of previousu calculated accleration
        for (int i =0; i <_N_real; i++){
            a_old[i][0] = a_new[i][0]; // when k = 0, a_new is the Newtownian term which calculated before
            a_old[i][1] = a_new[i][1];
            a_old[i][2] = a_new[i][2];
        }
        // now add on the non-constant term
        for (int i = 0; i < _N_real; i++){ // a_j is used to update a_i and vice versa
            double non_constx = 0.;
            double non_consty = 0.;
            double non_constz = 0.;
            for (int j = 0; j < _N_real; j++){
                if (j != i){
                    const double dxij = particles[i].x - particles[j].x;
                    const double dyij = particles[i].y - particles[j].y;
                    const double dzij = particles[i].z - particles[j].z;
                    const double r2ij = dxij*dxij + dyij*dyij + dzij*dzij;
                    const double rij = sqrt(r2ij);
                    non_constx += (G*particles[j].m*dxij/(r2ij*rij))*(dxij*a_old[j][0]+dyij*a_old[j][1]+dzij*a_old[j][2])/(2.*C*C)+\
                                        (7./(2.*C*C))*G*particles[j].m*a_old[j][0]/rij;
                    non_consty += (G*particles[j].m*dyij/(r2ij*rij))*(dxij*a_old[j][0]+dyij*a_old[j][1]+dzij*a_old[j][2])/(2.*C*C)+\
                                        (7./(2.*C*C))*G*particles[j].m*a_old[j][1]/rij;
                    non_constz += (G*particles[j].m*dzij/(r2ij*rij))*(dxij*a_old[j][0]+dyij*a_old[j][1]+dzij*a_old[j][2])/(2.*C*C)+\
                                        (7./(2.*C*C))*G*particles[j].m*a_old[j][2]/rij;
                }
            }
            a_new[i][0] = (a_const[i][0] + non_constx);
            a_new[i][1] = (a_const[i][1] + non_consty);
            a_new[i][2] = (a_const[i][2] + non_constz);
        }

        // break out loop if a_new is converging
        int breakout = 0;
        for (int i = 0; i < _N_real; i++){
            double dx = fabs(a_new[i][0] - a_old[i][0])/a_old[i][0];
            double dy = fabs(a_new[i][1] - a_old[i][1])/a_old[i][1];
            double dz = fabs(a_new[i][2] - a_old[i][2])/a_old[i][2];
            /*if (i ==1){FILE* d = fopen("difference.txt","a"); // this is used to check if the loop is giving resonable output
            fprintf(d, "number %d: %d: %.30e \n", k, i, dx);
            fclose(d);
            }*/
            if ((dx<1.e-30) && (dy <1.e-30) && (dz<1.e-30)){
                breakout += 1;
            }
        }
        if (breakout == 2){
            //printf(">>>>Precision reached in round %d<<<< \n", k);
            break;
        }
    }
    // update acceleration in particles
    for (int i = 0; i <_N_real;i++){
        particles[i].ax += a_new[i][0] - a_newton[i][0]; // substract newtonian term off since gravity routine already put it in
        particles[i].ay += a_new[i][1] - a_newton[i][1];
        particles[i].az += a_new[i][2] - a_newton[i][2];
    }
                    
}

double rebx_gr_full_hamiltonian(const struct reb_simulation* const sim, const struct rebx_params_gr_full* const params){
    const double C = params->c;
    const int _N_real = sim->N - sim->N_var;
    const double G = sim->G;
    struct reb_particle* const particles = sim->particles;

	double e_kin = 0.;
	double e_pot = 0.;
	double e_pn  = 0.;
	for (int i=0;i<_N_real;i++){
		struct reb_particle pi = particles[i];
		// Calculate p starts
        double pix=0.;
        double piy=0.;
        double piz=0.;

        double sumk = 0.;
        for (int k=0;k<_N_real;k++){
            if (k!=i){
                struct reb_particle pk = particles[k];
                double xik = pk.x - pi.x;
                double yik = pk.y - pi.y;
                double zik = pk.z - pi.z;
                double rik = sqrt(xik*xik + yik*yik + zik*zik);
                sumk -= 2.*G*pk.m/rik;
            }
        }
        
        double vi2 = pi.vx*pi.vx + pi.vy*pi.vy + pi.vz*pi.vz;
        double fac = pi.m*(1. + 0.5*vi2/C/C);

		for (int j=0;j<_N_real;j++){
			if (j != i){
				struct reb_particle pj = particles[j];
                double xij = pj.x - pi.x;
                double yij = pj.y - pi.y;
                double zij = pj.z - pi.z;
                double rij2 = xij*xij + yij*yij + zij*zij;
                double rij = sqrt(rij2);
                double rijdotvj = pj.vx*xij + pj.vy*yij + pj.vz*zij;
                double pfac = pi.m*pj.m/rij;

                pix += pfac*(6.*pi.vx - 7.*pj.vx - rijdotvj*xij/rij2);
                piy += pfac*(6.*pi.vy - 7.*pj.vy - rijdotvj*yij/rij2);
                piz += pfac*(6.*pi.vz - 7.*pj.vz - rijdotvj*zij/rij2);

                double rijdotvi = pi.vx*xij + pi.vy*yij + pi.vz*zij;
                double vidotvj = pi.vx*pj.vx + pi.vy*pj.vy + pi.vz*pj.vz;
                
                e_pn -= G/(4.*C*C)*pi.m*pj.m/rij*(6.*vi2 - 7*vidotvj - rijdotvi*rijdotvj/rij2 + sumk);
            }
        }

        pix *= G/(2.*C*C);
        piy *= G/(2.*C*C);
        piz *= G/(2.*C*C);
        
        pix += pi.vx*fac;
        piy += pi.vy*fac;
        piz += pi.vz*fac;
        // Calculate p ends
     
        double pi2 = pix*pix + piy*piy + piz*piz;
        e_kin += pi2/(2.*pi.m);
        e_pn -= pi.m/(8.*C*C)*vi2*vi2;

        for (int j=i+1; j<_N_real; j++){ // classic potential
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


