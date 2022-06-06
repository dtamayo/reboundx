/**
 * @file    stochastic_forces.c
 * @brief   Add stochastic forces
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE
 * Copyright (c) 2022 Hanno Rein, Dan Tamayo
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
 * $Stochastic Forces$     
 *
 * ======================= ===============================================
 * Authors                 H. Rein
 * Based on                `Rein and Papaloizou 2009 <https://ui.adsabs.harvard.edu/abs/2009A%26A...497..595R/abstract>`_.
 * Implementation Paper    `Rein and Choksi 2022 <https://iopscience.iop.org/article/10.3847/2515-5172/ac6e41>`_.
 * C Example               :ref:`c_example_stochastic_forces`
 * Python Example          `StochasticForces.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/StochasticForces.ipynb>`_, `StochasticForcesCartesian.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/StochasticForcesCartesian.ipynb>`_,
 * ======================= ===============================================
 * 
 * This applies stochastic forces to particles in the simulation.  
 * 
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 *
 * All particles which have the field kappa set, will experience stochastic forces.
 * The particle with index 0 cannot experience stochastic forces.
 *
 * ============================ =========== ==================================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================================
 * kappa (double)               Yes         Strength of stochastic forces relative to gravity from central object 
 * tau_kappa (double)           No          Auto-correlation time of stochastic forces. Defaults to orbital period if not set.
 *                                          The units are relative to the current orbital period.
 * ============================ =========== ==================================================================================
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "reboundx.h"


static void rebx_random_normal2(struct reb_simulation* r, double* n0, double* n1){
	double v1,v2,rsq=1.;
	while(rsq>=1. || rsq<1.0e-12){
		v1=2.*((double)rand_r(&(r->rand_seed)))/((double)(RAND_MAX))-1.0;
		v2=2.*((double)rand_r(&(r->rand_seed)))/((double)(RAND_MAX))-1.0;
		rsq=v1*v1+v2*v2;
	}
	*n0 = v1*sqrt(-2.*log(rsq)/rsq);
	*n1 = v2*sqrt(-2.*log(rsq)/rsq);
}


void rebx_stochastic_forces(struct reb_simulation* const sim, struct rebx_force* const radiation_forces, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    struct reb_particle com = particles[0];
    
    for (int i=0; i<N; i++){
        double* kappa = rebx_get_param(rebx, particles[i].ap, "kappa");
        if (i>0 && kappa != NULL){
            double* stochastic_force_r = rebx_get_param(rebx, particles[i].ap, "stochastic_force_r");
            if (stochastic_force_r == NULL) { // First run?
                rebx_set_param_double(rebx, (struct rebx_node**)&particles[i].ap, "stochastic_force_r", 0.);
                stochastic_force_r = rebx_get_param(rebx, particles[i].ap, "stochastic_force_r");
            }
            double* stochastic_force_phi = rebx_get_param(rebx, particles[i].ap, "stochastic_force_phi");
            if (stochastic_force_phi == NULL) { // First run?
                rebx_set_param_double(rebx, (struct rebx_node**)&particles[i].ap, "stochastic_force_phi", 0.);
                stochastic_force_phi = rebx_get_param(rebx, particles[i].ap, "stochastic_force_phi");
            }

            const struct reb_particle p = particles[i];

            // Get auto-correlation time
            int err=0;
            struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, particles[i], com, &err);
            if (err){
                reb_error(sim, "An error occured during the orbit calculation in rebx_stochastic_forces.\n");
                return;
            }
            double tau = o.P; // Default is current orbital period.
            
            double* tau_kappa = rebx_get_param(rebx, particles[i].ap, "tau_kappa");
            if (tau_kappa != NULL){
                tau *= *tau_kappa;
            }

            double dt = sim->dt_last_done;
            
            double prefac = exp(-dt/tau);

            // Decay
            *stochastic_force_r = (*stochastic_force_r) * prefac;
            *stochastic_force_phi = (*stochastic_force_phi) * prefac;

            double variance = 1.- prefac*prefac;
            if (variance <0.){
                reb_error(sim, "Timestep is larger than the correlation time for stochastic forces.\n");
                return;
            }
            double std = sqrt(variance);

            double n0, n1;
            rebx_random_normal2(sim, &n0, &n1);
            
            // Excitation
            *stochastic_force_r = (*stochastic_force_r) + n0*std;
            *stochastic_force_phi = (*stochastic_force_phi) + n1*std;

            const double dx = p.x - com.x; 
            const double dy = p.y - com.y;
            const double dz = p.z - com.z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz);
            
            const double dvx = p.vx - com.vx; 
            const double dvy = p.vy - com.vy;
            const double dvz = p.vz - com.vz;
            const double dv = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

            const double force_prefac = (*kappa) *sim->G/(dr*dr)*com.m;
            particles[i].ax += force_prefac*(*stochastic_force_r*dx/dr + *stochastic_force_phi*dvx/dv);
            particles[i].ay += force_prefac*(*stochastic_force_r*dy/dr + *stochastic_force_phi*dvy/dv);
            particles[i].az += force_prefac*(*stochastic_force_r*dz/dr + *stochastic_force_phi*dvz/dv);


		    com = reb_get_com_of_pair(com, p);
        }
        double* kappa_x = rebx_get_param(rebx, particles[i].ap, "kappa_x");
        if (kappa_x != NULL){
            double* stochastic_force_x = rebx_get_param(rebx, particles[i].ap, "stochastic_force_x");
            if (stochastic_force_x == NULL) { // First run?
                rebx_set_param_double(rebx, (struct rebx_node**)&particles[i].ap, "stochastic_force_x", 0.);
                stochastic_force_x = rebx_get_param(rebx, particles[i].ap, "stochastic_force_x");
            }
            
            double* tau_kappa_x = rebx_get_param(rebx, particles[i].ap, "tau_kappa_x");
            if (tau_kappa_x == NULL){
                reb_error(sim, "Need to set tau_kappa_x to enable stochastic forces.\n");
                return;
            }

            double dt = sim->dt_last_done;
            double prefac = exp(-dt/ (*tau_kappa_x));

            // Decay
            *stochastic_force_x = (*stochastic_force_x) * prefac;

            // Excitation
            double variance = 1.- prefac*prefac;
            if (variance <0.){
                reb_error(sim, "Timestep is larger than the correlation time for stochastic forces.\n");
                return;
            }
            double std = (*kappa_x)*sqrt(variance);
            double n0, n1;
            rebx_random_normal2(sim, &n0, &n1);
            *stochastic_force_x = (*stochastic_force_x) + n0*std;
            
            particles[i].ax += *stochastic_force_x;
        }
        double* kappa_y = rebx_get_param(rebx, particles[i].ap, "kappa_y");
        if (kappa_y != NULL){
            double* stochastic_force_y = rebx_get_param(rebx, particles[i].ap, "stochastic_force_y");
            if (stochastic_force_y == NULL) { // First run?
                rebx_set_param_double(rebx, (struct rebx_node**)&particles[i].ap, "stochastic_force_y", 0.);
                stochastic_force_y = rebx_get_param(rebx, particles[i].ap, "stochastic_force_y");
            }
            
            double* tau_kappa_y = rebx_get_param(rebx, particles[i].ap, "tau_kappa_y");
            if (tau_kappa_y == NULL){
                reb_error(sim, "Need to set tau_kappa_y to enable stochastic forces.\n");
                return;
            }

            double dt = sim->dt_last_done;
            double prefac = exp(-dt/ (*tau_kappa_y));

            // Decay
            *stochastic_force_y = (*stochastic_force_y) * prefac;

            // Excitation
            double variance = 1.- prefac*prefac;
            if (variance <0.){
                reb_error(sim, "Timestep is larger than the correlation time for stochastic forces.\n");
                return;
            }
            double std = (*kappa_y)*sqrt(variance);
            double n0, n1;
            rebx_random_normal2(sim, &n0, &n1);
            *stochastic_force_y = (*stochastic_force_y) + n0*std;
            
            particles[i].ay += *stochastic_force_y;
        }
        double* kappa_z = rebx_get_param(rebx, particles[i].ap, "kappa_z");
        if (kappa_z != NULL){
            double* stochastic_force_z = rebx_get_param(rebx, particles[i].ap, "stochastic_force_z");
            if (stochastic_force_z == NULL) { // First run?
                rebx_set_param_double(rebx, (struct rebx_node**)&particles[i].ap, "stochastic_force_z", 0.);
                stochastic_force_z = rebx_get_param(rebx, particles[i].ap, "stochastic_force_z");
            }
            
            double* tau_kappa_z = rebx_get_param(rebx, particles[i].ap, "tau_kappa_z");
            if (tau_kappa_z == NULL){
                reb_error(sim, "Need to set tau_kappa_z to enable stochastic forces.\n");
                return;
            }

            double dt = sim->dt_last_done;
            double prefac = exp(-dt/ (*tau_kappa_z));

            // Decay
            *stochastic_force_z = (*stochastic_force_z) * prefac;

            // Excitation
            double variance = 1.- prefac*prefac;
            if (variance <0.){
                reb_error(sim, "Timestep is larger than the correlation time for stochastic forces.\n");
                return;
            }
            double std = (*kappa_z)*sqrt(variance);
            double n0, n1;
            rebx_random_normal2(sim, &n0, &n1);
            *stochastic_force_z = (*stochastic_force_z) + n0*std;
            
            particles[i].az += *stochastic_force_z;
        }
    }
}

