/**
 * @file    gr.c
 * @brief   Post-newtonian general relativity corrections arising from a single massive body
 * @author  Pengshuai (Sam) Shi, Dan Tamayo, Hanno Rein <tamayo.daniel@gmail.com>
 *
 * @section     LICENSE
 * Copyright (c) 2015 Pengshuai (Sam) Shi, Dan Tamayo, Hanno Rein
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
 * Authors                 P. Shi, D. Tamayo, H. Rein
 * Implementation Paper    *In progress*
 * Based on                `Anderson et al. 1975 <http://labs.adsabs.harvard.edu/adsabs/abs/1975ApJ...200..221A/>`_.
 * C Example               :ref:`c_example_gr`
 * Python Example          `GeneralRelativity.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/GeneralRelativity.ipynb>`_.
 * ======================= ===============================================
 * 
 * This assumes that the masses are dominated by a single central body, and should be good enough for most applications with planets orbiting single stars.
 * It ignores terms that are smaller by of order the mass ratio with the central body.
 * It gets both the mean motion and precession correct, and will be significantly faster than :ref:`gr_full`, particularly with several bodies.
 * Adding this effect to several bodies is NOT equivalent to using gr_full.
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
 * If no particles have gr_source set, effect will assume the particle at index 0 in the particles array is the source.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * gr_source (int)              No          Flag identifying the particle as the source of perturbations.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_gr(struct reb_simulation* const sim, const double C2, const int source_index){
    const int N_real = sim->N - sim->N_var;
    const double G = sim->G;

    struct reb_particle* const ps = malloc(N_real*sizeof(*ps));
    struct reb_particle* const ps_j = malloc(N_real*sizeof(*ps_j));
    memcpy(ps, sim->particles, N_real*sizeof(*ps));
    
    // Calculate Newtonian accelerations 
    for(int i=0; i<N_real; i++){
        ps[i].ax = 0.;
        ps[i].ay = 0.;
        ps[i].az = 0.;
    }

    for(int i=0; i<N_real; i++){
        const struct reb_particle pi = ps[i];
        for(int j=i+1; j<N_real; j++){
            const struct reb_particle pj = ps[j];
            const double dx = pi.x - pj.x;
            const double dy = pi.y - pj.y;
            const double dz = pi.z - pj.z;
            const double softening2 = sim->softening*sim->softening;
            const double r2 = dx*dx + dy*dy + dz*dz + softening2;
            const double r = sqrt(r2);
            const double prefac = G/(r2*r);
            ps[i].ax -= prefac*pj.m*dx;
            ps[i].ay -= prefac*pj.m*dy;
            ps[i].az -= prefac*pj.m*dz;
            ps[j].ax += prefac*pi.m*dx;
            ps[j].ay += prefac*pi.m*dy;
            ps[j].az += prefac*pi.m*dz;
        }
    }
   
    // Transform to Jacobi coordinates
    const struct reb_particle source = ps[0];
	const double mu = G*source.m;
    double* const eta = malloc(N_real*sizeof(*eta));
    eta[0] = ps[0].m;
    for (unsigned int i=1;i<N_real;i++){
        eta[i] = eta[i-1] + ps[i].m;
    }
  
    to_jacobi_posvel(ps, ps_j, eta, ps, N_real);
    to_jacobi_acc(ps, ps_j, eta, ps, N_real);

    for (int i=1; i<N_real; i++){
        struct reb_particle p = ps_j[i];
        struct reb_vec3d vi;
        vi.x = p.vx;
        vi.y = p.vy;
        vi.z = p.vz;
        double vi2, A;
        const double ri = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        for(int q=0; q<10; q++){
            vi2 = vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
            A = (0.5*vi2 + 3.*mu/ri)/C2;
            vi.x = p.vx/(1.-A);
            vi.y = p.vy/(1.-A);
            vi.z = p.vz/(1.-A);
        }

        const double B = (mu/ri - 1.5*vi2)*mu/(ri*ri*ri)/C2;
        const double vdota = vi.x*p.ax + vi.y*p.ay + vi.z*p.az;
        const double vdotr = vi.x*p.x + vi.y*p.x + vi.z*p.z;
        const double rdotrdot = p.x*p.vx + p.y*p.vy + p.z*p.vz;
        const double D = (vdota + B*vdotr - 3.*mu/(ri*ri*ri)*rdotrdot)/C2;

        ps_j[i].ax = B*(1.-A)*p.x - A*p.ax - D*vi.x;
        ps_j[i].ay = B*(1.-A)*p.y - A*p.ay - D*vi.y;
        ps_j[i].az = B*(1.-A)*p.z - A*p.az - D*vi.z;
    }
    ps_j[0].ax = 0.;
    ps_j[0].ay = 0.;
    ps_j[0].az = 0.;

    to_inertial_acc(ps, ps_j, eta, ps, N_real);
    for (int i=0; i<N_real; i++){
        sim->particles[i].ax += ps[i].ax;
        sim->particles[i].ay += ps[i].ay;
        sim->particles[i].az += ps[i].az;
    }
}

/*static void rebx_calculate_gr(struct reb_simulation* const sim, const double C2, const int source_index){
    const int N_real = sim->N - sim->N_var;
    const double G = sim->G;
    struct reb_particle* const ps = sim->particles;
    

    // Calculate Newtonian accelerations 
    struct reb_vec3d* const a = calloc(N_real, sizeof(*a));
    for(int i=0; i<N_real; i++){
        const struct reb_particle pi = ps[i];
        for(int j=i+1; j<N_real; j++){
            const struct reb_particle pj = ps[j];
            const double dx = pi.x - pj.x;
            const double dy = pi.y - pj.y;
            const double dz = pi.z - pj.z;
            const double softening2 = sim->softening*sim->softening;
            const double r2 = dx*dx + dy*dy + dz*dz + softening2;
            const double r = sqrt(r2);
            const double prefac = G/(r2*r);
            a[i].x -= prefac*pj.m*dx;
            a[i].y -= prefac*pj.m*dy;
            a[i].z -= prefac*pj.m*dz;
            a[j].x += prefac*pi.m*dx;
            a[j].y += prefac*pi.m*dy;
            a[j].z += prefac*pi.m*dz;
        }
    }
    
    const struct reb_particle source = ps[source_index];
	const double mu = G*source.m;
   
    // calculate pseudovelocities
    double* const const_term = malloc(N_real*sizeof(*const_term));
    double Mtot = 0.;
    for (int i=0; i<N_real; i++){
        Mtot += ps[i].m;
        if(i == source_index){
            continue;
        }
		const struct reb_particle pi = ps[i];
		const double dx = pi.x - source.x;
		const double dy = pi.y - source.y;
		const double dz = pi.z - source.z;
		const double r2 = dx*dx + dy*dy + dz*dz;
		const_term[i] = 3.*mu/sqrt(r2);
    }		
   
    struct reb_vec3d* const v_tilde = malloc(N_real*sizeof(*v_tilde));
    struct reb_vec3d* const nu_tilde = malloc(N_real*sizeof(*nu_tilde));
    for(int i=0; i<N_real; i++){
        v_tilde[i].x = ps[i].vx;
        v_tilde[i].y = ps[i].vy;
        v_tilde[i].z = ps[i].vz;
    }
    
    double v_tilde1x = ps[1].vx/(1.-(ps[1].vx*ps[1].vx/2. + const_term[1])/C2);
    fprintf(stderr, "%e\n***\n", ps[1].vx - v_tilde1x);

    for(int q=0; q<10; q++){
        // Calculate new v_tilde for COM
        struct reb_vec3d v_tildeCOM = {0.};
        for(int i=0; i<N_real; i++){
            v_tildeCOM.x += ps[i].m*v_tilde[i].x;
            v_tildeCOM.y += ps[i].m*v_tilde[i].y;
            v_tildeCOM.z += ps[i].m*v_tilde[i].z;
        }
        v_tildeCOM.x /= Mtot;
        v_tildeCOM.y /= Mtot;
        v_tildeCOM.z /= Mtot;
        fprintf(stderr, "v_tildeCOM.x : %e\n", v_tildeCOM.x);
        // Calculate new nu = v_tilde - v_tildeCOM
        for(int i=0; i<N_real; i++){
            nu_tilde[i].x = v_tilde[i].x - v_tildeCOM.x;
            nu_tilde[i].y = v_tilde[i].y - v_tildeCOM.y;
            nu_tilde[i].z = v_tilde[i].z - v_tildeCOM.z;
        } 
        
        struct reb_vec3d v_PNCOM = {0.};
        for (int i=0; i<N_real; i++){
            if(i == source_index){
                continue;
            }
            const struct reb_particle pi = ps[i];
            const double nu_tilde2 = nu_tilde[i].x*nu_tilde[i].x + nu_tilde[i].y*nu_tilde[i].y + nu_tilde[i].z*nu_tilde[i].z;
            const double prefac = (nu_tilde2/2. + const_term[i])/C2;
            const double vPNx = -prefac*nu_tilde[i].x;
            const double vPNy = -prefac*nu_tilde[i].y;
            const double vPNz = -prefac*nu_tilde[i].z;
            v_PNCOM.x += vPNx;
            v_PNCOM.y += vPNy;
            v_PNCOM.x += vPNz;
            v_tilde[i].x = ps[i].vx - vPNx;
            v_tilde[i].y = ps[i].vy - vPNy;
            v_tilde[i].z = ps[i].vz - vPNz;
        }		
        v_PNCOM.x /= Mtot;
        v_PNCOM.y /= Mtot;
        v_PNCOM.z /= Mtot;
        for(int i=0; i<N_real; i++){
            v_tilde[i].x += ps[i].m*v_PNCOM.x;
            v_tilde[i].y += ps[i].m*v_PNCOM.y;
            v_tilde[i].z += ps[i].m*v_PNCOM.z;
        }
        fprintf(stderr, "%d\t%e\n", q, ps[1].vx - v_tilde[1].x);
        
    }
*/
    /*
        const double vx = pi.vx;
		const double vy = pi.vy;
		const double vz = pi.vz;
		const double v2 = vx*vx + vy*vy + vz*vz;
        
        const double a1_x = (mu*mu*dx/(r2*r2) - 3.*mu*v2*dx/(2.*r2*r))/C2;
		const double a1_y = (mu*mu*dy/(r2*r2) - 3.*mu*v2*dy/(2.*r2*r))/C2;
		const double a1_z = (mu*mu*dz/(r2*r2) - 3.*mu*v2*dz/(2.*r2*r))/C2;

		const double va = vx*ax + vy*ay + vz*az;
		const double rv = dx*vx + dy*vy + dz*vz;
	
        ax = a1_x-(va*vx + v2*ax/2. + 3.*mu*(ax*r-vx*rv/r)/r2)/C2;
        ay = a1_y-(va*vy + v2*ay/2. + 3.*mu*(ay*r-vy*rv/r)/r2)/C2;
        az = a1_z-(va*vz + v2*az/2. + 3.*mu*(az*r-vz*rv/r)/r2)/C2;

		particles[i].ax += ax;
		particles[i].ay += ay; 
		particles[i].az += az;
		particles[source_index].ax -= pi.m/source.m*ax; 
		particles[source_index].ay -= pi.m/source.m*ay; 
		particles[source_index].az -= pi.m/source.m*az; 
    }

}
*/
/*
static void rebx_calculate_gr(struct reb_simulation* const sim, const double C2, const int source_index){
    const int N_real = sim->N - sim->N_var;
    const double G = sim->G;
    struct reb_particle* const particles = sim->particles;
    const unsigned int _gravity_ignore_10 = sim->gravity_ignore_terms==1;
    
	const double mu = G*particles[source_index].m;
    double aoverm10x, aoverm10y, aoverm10z;

    if (_gravity_ignore_10){
        const double dx = particles[0].x - particles[1].x;
        const double dy = particles[0].y - particles[1].y;
        const double dz = particles[0].z - particles[1].z;
        const double softening2 = sim->softening*sim->softening;
        const double r2 = dx*dx + dy*dy + dz*dz + softening2;
        const double r = sqrt(r2);
        const double prefac = G/(r2*r);
        
        aoverm10x = prefac*dx;
        aoverm10y = prefac*dy;
        aoverm10z = prefac*dz;
    }
	
    const struct reb_particle source = particles[source_index];
    for (int i=0; i<N_real; i++){
        if(i == source_index){
            continue;
        }
		const struct reb_particle pi = particles[i];
		
		const double dx = pi.x - source.x;
		const double dy = pi.y - source.y;
		const double dz = pi.z - source.z;
		const double r2 = dx*dx + dy*dy + dz*dz;
		const double r = sqrt(r2);
		const double vx = pi.vx;
		const double vy = pi.vy;
		const double vz = pi.vz;
		const double v2 = vx*vx + vy*vy + vz*vz;
        double ax = pi.ax;
        double ay = pi.ay;
        double az = pi.az;
        if(_gravity_ignore_10 && i==1){
            ax += particles[0].m*aoverm10x;
            ay += particles[0].m*aoverm10y;
            az += particles[0].m*aoverm10z;
        }
        if(_gravity_ignore_10 && i==0){
            ax -= particles[1].m*aoverm10x;
            ay -= particles[1].m*aoverm10y;
            az -= particles[1].m*aoverm10z;
        }
        
        const double a1_x = (mu*mu*dx/(r2*r2) - 3.*mu*v2*dx/(2.*r2*r))/C2;
		const double a1_y = (mu*mu*dy/(r2*r2) - 3.*mu*v2*dy/(2.*r2*r))/C2;
		const double a1_z = (mu*mu*dz/(r2*r2) - 3.*mu*v2*dz/(2.*r2*r))/C2;

		const double va = vx*ax + vy*ay + vz*az;
		const double rv = dx*vx + dy*vy + dz*vz;
	
        ax = a1_x-(va*vx + v2*ax/2. + 3.*mu*(ax*r-vx*rv/r)/r2)/C2;
        ay = a1_y-(va*vy + v2*ay/2. + 3.*mu*(ay*r-vy*rv/r)/r2)/C2;
        az = a1_z-(va*vz + v2*az/2. + 3.*mu*(az*r-vz*rv/r)/r2)/C2;

		particles[i].ax += ax;
		particles[i].ay += ay; 
		particles[i].az += az;
		particles[source_index].ax -= pi.m/source.m*ax; 
		particles[source_index].ay -= pi.m/source.m*ay; 
		particles[source_index].az -= pi.m/source.m*az; 
    }	
}
*/
void rebx_gr(struct reb_simulation* const sim, struct rebx_effect* const gr){ // First find gr sources
    double* c = rebx_get_param_check(gr, "c", REBX_TYPE_DOUBLE);
    if (c == NULL){
        reb_error(sim, "Need to set speed of light in gr effect.  See examples in documentation.\n");
        return;
    }
    const double C2 = (*c)*(*c);
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    for (int i=0; i<N_real; i++){
        if (rebx_get_param_check(&particles[i], "gr_source", REBX_TYPE_INT) != NULL){
            rebx_calculate_gr(sim, C2, i);
            return;                 // only apply effect for first gr_source found.  For multiple sources, need gr_full
        }
    }
    rebx_calculate_gr(sim, C2, 0);  // gr_source not found, default to index=0
}

static double rebx_calculate_gr_hamiltonian(struct reb_simulation* const sim, const double C2, const int source_index){
    const int N_real = sim->N - sim->N_var;
    const double G = sim->G;

    struct reb_particle* const ps_j = malloc(N_real*sizeof(*ps_j));
    struct reb_particle* const ps = sim->particles; 
    // Calculate Newtonian potentials

    double V_newt = 0.;
    for(int i=0; i<N_real; i++){
        const struct reb_particle pi = ps[i];
        for(int j=i+1; j<N_real; j++){
            const struct reb_particle pj = ps[j];
            const double dx = pi.x - pj.x;
            const double dy = pi.y - pj.y;
            const double dz = pi.z - pj.z;
            const double softening2 = sim->softening*sim->softening;
            const double r2 = dx*dx + dy*dy + dz*dz + softening2;
            V_newt -= G*pi.m*pj.m/sqrt(r2);
        }
    }
   
    // Transform to Jacobi coordinates
    const struct reb_particle source = ps[0];
	const double mu = G*source.m;
    double* const eta = malloc(N_real*sizeof(*eta));
    double* const m_j = malloc(N_real*sizeof(*m_j));
    eta[0] = ps[0].m;
    for (unsigned int i=1;i<N_real;i++){
        eta[i] = eta[i-1] + ps[i].m;
        m_j[i] = ps[i].m*eta[i-1]/eta[i];
    }
    m_j[0] = eta[N_real-1];

    to_jacobi_posvel(ps, ps_j, eta, ps, N_real);

    double T = 0.5*m_j[0]*(ps_j[0].vx*ps_j[0].vx + ps_j[0].vy*ps_j[0].vy + ps_j[0].vz*ps_j[0].vz);
    double V_PN = 0.;
    for (int i=1; i<N_real; i++){
        struct reb_particle p = ps_j[i];
        const double rdoti2 = p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
        double vi2 = rdoti2;
        double A;
        const double ri = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        for(int q=0; q<10; q++){
            A = (0.5*vi2 + 3.*mu/ri)/C2;
            vi2 = rdoti2/((1.-A)*(1.-A));
        }

        V_PN += m_j[i]*(0.5*mu*mu/(ri*ri) - 0.125*vi2*vi2 - 1.5*mu*vi2/ri);
        T += 0.5*m_j[i]*vi2;
    }
    V_PN /= C2;
    
	return T + V_newt + V_PN;
}

double rebx_gr_hamiltonian(struct reb_simulation* const sim, const struct rebx_effect* const gr){ 
    double* c = rebx_get_param_check(gr, "c", REBX_TYPE_DOUBLE);
    if (c == NULL){
        reb_error(sim, "Need to set speed of light in gr effect.  See examples in documentation.\n");
        return 0;
    }
    const double C2 = (*c)*(*c);
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    for (int i=0; i<N_real; i++){
        if (rebx_get_param_check(&particles[i], "gr_source", REBX_TYPE_INT) != NULL){
            return rebx_calculate_gr_hamiltonian(sim, C2, i);
        }
    }
    return rebx_calculate_gr_hamiltonian(sim, C2, 0);
}

