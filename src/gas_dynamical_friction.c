/**
 * @file    gas_dynamical_friction.c
 * @brief   Gas drag from a thin, disk with a power-law density profile 
 * @author  Aleksey Generozov
 * 
 * 
 *
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
 * $Gas Effects$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 A. Generozov, H. Perets
 * Implementation Paper    `Generozov and Perets 2022 <https://arxiv.org/abs/2212.11301>`_
 * Based on                `Ostriker 1999 (with simplifications) <https://ui.adsabs.harvard.edu/abs/1999ApJ...513..252O/abstract>`_, `Just et al 2012 <https://ui.adsabs.harvard.edu/abs/2012ApJ...758...51J/abstract>`_.
 * C Example               :ref:`c_example_gas_dynamical_friction`
 * Python Example          `gas_dynamical_friction.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/gas_dynamical_friction.ipynb>`_
 *                        
 * 
 * ======================= ===============================================
 * 
 * 
 * **Effect Parameters**
 * 
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * rhog (double)                Yes         Normalization of density. Density in the disk midplane is rhog*r^alpha_rhog
 * alpha_rhog (double)          Yes         Power-law slope of the power-law density profile.
 * cs (double)                  Yes         Normalization of the sound speed. Sound speed has profile cs*r^alpha_cs
 * alpha_cs (double)            Yes         Power-law slope of the sound speed
 * xmin (double)                Yes         Dimensionless parameter that determines the Coulomb logarithm (ln(L) =log (1/xmin))
 * hr (double)                  Yes         Aspect ratio of the disk
 * Qd (double)                  Yes         Prefactor for geometric drag
 * ============================ =========== ==================================================================
 * 
 * 
 * **Particle Parameters**
 * 
 * None.
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "reboundx.h"

static double mach_piece_sub(const double mach){
    //Using powerlaw expansion at low mach numbers for numerical reasons
    if (mach<0.02){
        return mach*mach*mach/3.+mach*mach*mach*mach*mach/5.;
    }
    //Subsonic expression from Ostriker...
    return 0.5*log((1.0+mach)/(1.0-mach))-mach;

}


static double calculate_pre_factor(const double mach, const double vel, const double t,
    const double xmin){
    //Simplified version of the Ostriker dynamical friction formula...
    const double coul=log(1.0/xmin);
    if (mach>=1.0){
        return coul;
    }
    else{
        return fmin(coul, mach_piece_sub(mach));
    }

}

static void get_vrel_disk(const struct reb_particle p, const double GMBH, double* vrel, const double hr){
    const double rcyl = sqrt(p.x*p.x+p.y*p.y);
    const double sin_phi= p.y/rcyl;
    const double cos_phi= p.x/rcyl;

    const double vk=sqrt(GMBH/rcyl)*(1.0-hr*hr);
    vrel[0]=p.vx+vk*sin_phi;
    vrel[1]=p.vy-vk*cos_phi;
    vrel[2]=p.vz;

}

static void rebx_calculate_gas_dynamical_friction(struct reb_simulation* const sim, struct reb_particle* const particles,\
    const int N, const double rhog, const double alpha_rhog, const double cs, const double alpha_cs, const double xmin, const double hr, const double Qd){

    const int _N_real = sim->N - sim->N_var;
    const struct reb_particle bh = particles[0];
#pragma omp parallel for
    for (int i=1;i<_N_real;i++){
        const struct reb_particle p = particles[i];
        struct reb_particle diff = p;
        struct reb_particle bh2  = bh;
        reb_particle_isub(&diff, &bh2);
        const double rcyl=sqrt(diff.x*diff.x+diff.y*diff.y);
        double vrel [3];
        get_vrel_disk(diff, sim->G*bh.m, vrel, hr);
        const double vrel_norm=sqrt(vrel[0]*vrel[0]+vrel[1]*vrel[1]+vrel[2]*vrel[2]);
        const double mach=vrel_norm/(cs*pow(rcyl, alpha_cs));
        const double t=sim->t;

        const double integ=calculate_pre_factor(mach, vrel_norm, t, xmin);
        //Accounting for vertical dependence of the density with a Gaussian function
        //scale height is defined by user-defined aspect ratio. Truncate the disc vertically
        //at 10 scale heights...
        const double h=hr*rcyl;
        const double vert=(abs(diff.z)<(10*h))?exp(-diff.z*diff.z/(2.0*h*h)):0;
        const double rhog_loc=rhog*pow(rcyl, alpha_rhog)*vert;
        const double mp = p.m;
        const double rstar = p.r;
        const double fc=4.*M_PI*(sim->G*sim->G)*mp*(rhog_loc)/(vrel_norm*vrel_norm*vrel_norm)*integ+\
M_PI*rhog_loc*rstar*rstar*vrel_norm*Qd/mp;

        particles[i].ax -= fc*vrel[0];
        particles[i].ay -= fc*vrel[1];
        particles[i].az -= fc*vrel[2];

    }
}

void rebx_gas_dynamical_friction(struct reb_simulation* const sim, struct rebx_force* const force,\
    struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    double* rhog= rebx_get_param(rebx, force->ap, "gas_df_rhog");
    if (rhog == NULL){
        reb_simulation_error(sim, "Need to specify a gas density\n");
    }
    double* alpha_rhog= rebx_get_param(rebx, force->ap, "gas_df_alpha_rhog");
    if (alpha_rhog== NULL){
        reb_simulation_error(sim, "Need to specify a profile for gas density\n");
    }
    double* cs= rebx_get_param(rebx, force->ap, "gas_df_cs");
    if (cs == NULL){
        reb_simulation_error(sim, "Need to set a sound speed.\n");
    }
    double* alpha_cs= rebx_get_param(rebx, force->ap, "gas_df_alpha_cs");
    if (alpha_cs== NULL){
        reb_simulation_error(sim, "Need to specify a profile for the sound speed\n");
    }
    double* xmin= rebx_get_param(rebx, force->ap, "gas_df_xmin");
    if (xmin == NULL){
        reb_simulation_error(sim, "Need to set a cutoff.\n");
    }
    double* hr= rebx_get_param(rebx, force->ap, "gas_df_hr");
    if (hr == NULL){
        reb_simulation_error(sim, "Need an aspect ratio.\n");
    }
    double* Qd=rebx_get_param(rebx, force->ap, "gas_df_Qd");
    if (Qd == NULL){
        reb_simulation_error(sim, "Need to specify Qd");
    }

    rebx_calculate_gas_dynamical_friction(sim, particles, N, *rhog, *alpha_rhog, *cs, *alpha_cs, *xmin, *hr, *Qd);

}


