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
 * c (double)                   Yes         Speed of light, needs to be specified in the units used for the simulation.
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

static void reb_particles_transform_inertial_to_barycentric_posvelacc(const struct reb_particle* const particles, struct reb_particle* const p_b, const unsigned int N, const unsigned int N_active){
    double x0  = 0.;
    double y0  = 0.;
    double z0  = 0.;
    double vx0 = 0.;
    double vy0 = 0.;
    double vz0 = 0.;
    double ax0 = 0.;
    double ay0 = 0.;
    double az0 = 0.;
    double m0  = 0.;
#pragma omp parallel for reduction(+:x0) reduction(+:y0) reduction(+:z0) reduction(+:vx0) reduction(+:vy0) reduction(+:vz0) reduction(+:ax0) reduction(+:ay0) reduction(+:az0) reduction(+:m0)
    for (unsigned int i=0;i<N_active;i++){
        double m = particles[i].m;
        x0  += particles[i].x *m;
        y0  += particles[i].y *m;
        z0  += particles[i].z *m;
        vx0 += particles[i].vx*m;
        vy0 += particles[i].vy*m;
        vz0 += particles[i].vz*m;
        ax0 += particles[i].ax*m;
        ay0 += particles[i].ay*m;
        az0 += particles[i].az*m;
        m0  += m;
    }
    p_b[0].x  = x0/m0;
    p_b[0].y  = y0/m0;
    p_b[0].z  = z0/m0;
    p_b[0].vx = vx0/m0;
    p_b[0].vy = vy0/m0;
    p_b[0].vz = vz0/m0;
    p_b[0].ax = ax0/m0;
    p_b[0].ay = ay0/m0;
    p_b[0].az = az0/m0;
    p_b[0].m = m0;
    
#pragma omp parallel for 
    for (unsigned int i=1;i<N;i++){
        p_b[i].x  = particles[i].x  - p_b[0].x;
        p_b[i].y  = particles[i].y  - p_b[0].y;
        p_b[i].z  = particles[i].z  - p_b[0].z;
        p_b[i].vx = particles[i].vx - p_b[0].vx;
        p_b[i].vy = particles[i].vy - p_b[0].vy;
        p_b[i].vz = particles[i].vz - p_b[0].vz;
        p_b[i].ax = particles[i].ax - p_b[0].ax;
        p_b[i].ay = particles[i].ay - p_b[0].ay;
        p_b[i].az = particles[i].az - p_b[0].az;
        p_b[i].m  = particles[i].m;
    }
}

static void reb_particles_transform_barycentric_to_inertial_acc(struct reb_particle* const particles, const struct reb_particle* const p_b, const unsigned int N, const unsigned int N_active){
    const double mtot = p_b[0].m;
    double ax0  = 0.;
    double ay0  = 0.;
    double az0  = 0.;
#pragma omp parallel for reduction(+:ax0) reduction(+:ay0) reduction(+:az0) 
    for (unsigned int i=1;i<N_active;i++){
        double m = p_b[i].m;
        ax0 += p_b[i].ax*m/mtot;
        ay0 += p_b[i].ay*m/mtot;
        az0 += p_b[i].az*m/mtot;
        particles[i].m = m; // in case of merger/mass change
    }
    particles[0].ax  = p_b[0].ax - ax0;
    particles[0].ay  = p_b[0].ay - ay0;
    particles[0].az  = p_b[0].az - az0;
#pragma omp parallel for 
    for (unsigned int i=1;i<N;i++){
        particles[i].ax = p_b[i].ax+particles[0].ax;
        particles[i].ay = p_b[i].ay+particles[0].ay;
        particles[i].az = p_b[i].az+particles[0].az;
    }
}

static void rebx_calculate_gr_full(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double C2, const double G, const int max_iterations, const int gravity_ignore_10){
    
    double a_const[N][3]; // array that stores the value of the constant term
    struct reb_particle* const ps_b = malloc(N*sizeof(*ps_b));
    memcpy(ps_b, particles, N*sizeof(*ps_b));

    // Calculate Newtonian accelerations 
    for(int i=0; i<N; i++){
        ps_b[i].ax = 0.;
        ps_b[i].ay = 0.;
        ps_b[i].az = 0.;
    }

    for(int i=0; i<N; i++){
        const struct reb_particle pi = ps_b[i];
        for(int j=i+1; j<N; j++){
            const struct reb_particle pj = ps_b[j];
            const double dx = pi.x - pj.x;
            const double dy = pi.y - pj.y;
            const double dz = pi.z - pj.z;
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double r = sqrt(r2);
            const double prefac = G/(r2*r);
            ps_b[i].ax -= prefac*pj.m*dx;
            ps_b[i].ay -= prefac*pj.m*dy;
            ps_b[i].az -= prefac*pj.m*dz;
            ps_b[j].ax += prefac*pi.m*dx;
            ps_b[j].ay += prefac*pi.m*dy;
            ps_b[j].az += prefac*pi.m*dz;
        }
    }

    // Transform to barycentric coordinates
    struct reb_particle com = reb_simulation_com(sim);
    for (int i=0; i<N; i++){
        reb_particle_isub(&ps_b[i], &com);
    }
    for (int i=0; i<N; i++){
        // then compute the constant terms:
        double a_constx = 0.;
        double a_consty = 0.;
        double a_constz = 0.;
        // 1st constant part
        for (int j = 0; j< N; j++){
            if (j != i){
                const double dxij = ps_b[i].x - ps_b[j].x;
                const double dyij = ps_b[i].y - ps_b[j].y;
                const double dzij = ps_b[i].z - ps_b[j].z;
                const double rij2 = dxij*dxij + dyij*dyij + dzij*dzij;
                const double rij = sqrt(rij2);
                const double rij3 = rij2*rij;
                
                double a1 = 0.;
                for (int k = 0; k< N; k++){
                    if (k != i){
                        const double dxik = ps_b[i].x - ps_b[k].x;
                        const double dyik = ps_b[i].y - ps_b[k].y;
                        const double dzik = ps_b[i].z - ps_b[k].z;
                        const double rik = sqrt(dxik*dxik + dyik*dyik + dzik*dzik);
                        a1 += (4./(C2)) * G*particles[k].m/rik;
                    }
                }

                double a2 = 0.;
                for (int l = 0; l< N; l++){
                    if (l != j){
                        const double dxlj = ps_b[l].x - ps_b[j].x;
                        const double dylj = ps_b[l].y - ps_b[j].y;
                        const double dzlj = ps_b[l].z - ps_b[j].z;
                        const double rlj = sqrt(dxlj*dxlj + dylj*dylj + dzlj*dzlj);
                        a2 += (1./(C2)) * G*particles[l].m/rlj;
                    }
                }

                double a3;
                double vi2 = ps_b[i].vx*ps_b[i].vx + ps_b[i].vy*ps_b[i].vy + ps_b[i].vz*ps_b[i].vz;
                a3 = -vi2/(C2);

                double a4;
                double vj2 = ps_b[j].vx*ps_b[j].vx + ps_b[j].vy*ps_b[j].vy + ps_b[j].vz*ps_b[j].vz;
                a4 = -2.*vj2/(C2);

                double a5;
                a5 = (4./(C2)) * (ps_b[i].vx*ps_b[j].vx + ps_b[i].vy*ps_b[j].vy + ps_b[i].vz*ps_b[j].vz); 
                
                double a6;
                double a6_0 = dxij*ps_b[j].vx + dyij*ps_b[j].vy + dzij*ps_b[j].vz;
                a6 = (3./(2.*C2)) * a6_0*a6_0/rij2;
               
                double a7; // Newtonian piece of first ddot(r) piece
                a7 = (dxij*ps_b[j].ax+dyij*ps_b[j].ay+dzij*ps_b[j].az)/(2.*C2);
                
                double factor1 = a1 + a2 + a3 + a4 + a5 + a6 + a7;
                 
                a_constx += G*particles[j].m*dxij*factor1/rij3;
                a_consty += G*particles[j].m*dyij*factor1/rij3;
                a_constz += G*particles[j].m*dzij*factor1/rij3;
        
                // 2nd constant part
                
                const double dvxij = ps_b[i].vx - ps_b[j].vx;
                const double dvyij = ps_b[i].vy - ps_b[j].vy;
                const double dvzij = ps_b[i].vz - ps_b[j].vz;
                    
                double factor2 = dxij*(4.*ps_b[i].vx-3.*ps_b[j].vx)+dyij*(4.*ps_b[i].vy-3.*ps_b[j].vy)+dzij*(4.*ps_b[i].vz-3.*ps_b[j].vz);

                a_constx += G*particles[j].m/C2*(factor2*dvxij/rij3 + 7./2.*ps_b[j].ax/rij);
                a_consty += G*particles[j].m/C2*(factor2*dvyij/rij3 + 7./2.*ps_b[j].ay/rij);
                a_constz += G*particles[j].m/C2*(factor2*dvzij/rij3 + 7./2.*ps_b[j].az/rij);
            }
        }  

        a_const[i][0] = a_constx;
        a_const[i][1] = a_consty;
        a_const[i][2] = a_constz;
    }
    for (int i = 0; i <N; i++){
        ps_b[i].ax = a_const[i][0];
        ps_b[i].ay = a_const[i][1];
        ps_b[i].az = a_const[i][2];
    }

    // Now running the substitution again and again through the loop below
    for (int k=0; k<10; k++){ // you can set k as how many substitution you want to make
        double a_old[N][3]; // initialize an arry that stores the information of previousu calculated accleration
        for (int i =0; i <N; i++){
            a_old[i][0] = ps_b[i].ax;
            a_old[i][1] = ps_b[i].ay;
            a_old[i][2] = ps_b[i].az;
        }
        // now add on the non-constant term
        for (int i = 0; i < N; i++){ // a_j is used to update a_i and vice versa
            double non_constx = 0.;
            double non_consty = 0.;
            double non_constz = 0.;
            for (int j = 0; j < N; j++){
                if (j != i){
                    const double dxij = ps_b[i].x - ps_b[j].x;
                    const double dyij = ps_b[i].y - ps_b[j].y;
                    const double dzij = ps_b[i].z - ps_b[j].z;
                    const double rij = sqrt(dxij*dxij + dyij*dyij + dzij*dzij);
                    const double rij3 = rij*rij*rij;
                    const dotproduct = dxij*ps_b[j].ax+dyij*ps_b[j].ay+dzij*ps_b[j].az;

                    non_constx += (G*particles[j].m*dxij/rij3)*dotproduct/(2.*C2) + (7./(2.*C2))*G*particles[j].m*ps_b[j].ax/rij;
                    non_consty += (G*particles[j].m*dyij/rij3)*dotproduct/(2.*C2) + (7./(2.*C2))*G*particles[j].m*ps_b[j].ay/rij;
                    non_constz += (G*particles[j].m*dzij/rij3)*dotproduct/(2.*C2) + (7./(2.*C2))*G*particles[j].m*ps_b[j].az/rij;
                }
            }
            ps_b[i].ax = a_const[i][0] + non_constx;
            ps_b[i].ay = a_const[i][1] + non_consty;
            ps_b[i].az = a_const[i][2] + non_constz;
        }
        
        // break out loop if ps_b is converging
        double maxdev = 0.;
        double dx, dy, dz;
        for (int i = 0; i < N; i++){
            dx = (fabs(ps_b[i].ax) < DBL_EPSILON) ? 0. : fabs((ps_b[i].ax - a_old[i][0])/ps_b[i].ax);
            dy = (fabs(ps_b[i].ay) < DBL_EPSILON) ? 0. : fabs((ps_b[i].ay - a_old[i][1])/ps_b[i].ay);
            dz = (fabs(ps_b[i].az) < DBL_EPSILON) ? 0. : fabs((ps_b[i].az - a_old[i][2])/ps_b[i].az);
            
            if (dx > maxdev) { maxdev = dx; }
            if (dy > maxdev) { maxdev = dy; }
            if (dz > maxdev) { maxdev = dz; }
        }
        
        if (maxdev < DBL_EPSILON){
            break;
        }
        if (k==9){
            reb_simulation_warning(sim, "10 loops in rebx_gr_full did not converge.\n");
            fprintf(stderr, "Fractional Error: %e\n", maxdev);
        }
    }
   
    for (int i=0; i<N; i++){
        particles[i].ax += ps_b[i].ax;
        particles[i].ay += ps_b[i].ay;
        particles[i].az += ps_b[i].az;
    }
    
    free(ps_b);

}

void rebx_gr_full(struct reb_simulation* const sim, struct rebx_force* const gr_full, struct reb_particle* const particles, const int N){
    double* c = rebx_get_param(sim->extras, gr_full->ap, "c");
    if (c == NULL){
        reb_simulation_error(sim, "REBOUNDx Error: Need to set speed of light in gr effect.  See examples in documentation.\n");
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
        reb_simulation_error(sim, "Need to set speed of light in gr effect.  See examples in documentation.\n");
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
