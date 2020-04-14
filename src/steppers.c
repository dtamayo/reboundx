/**
 * @file    steppers.c
 * @brief   Wrappers to run REBOUND steps with different integrators
 * @author  Dan Tamayo, Hanno Rein <tamayo.daniel@gmail.com>
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
 * $Integration Steppers$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo, H. Rein
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
 * Based on                `Rein and Liu, 2012 <https://ui.adsabs.harvard.edu/abs/2012A%26A...537A.128R/abstract>`_.
 * C Example               None
 * Python Example          `CustomSplittingIntegrationSchemes.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/CustomSplittingIntegrationSchemes.ipynb>`_.
 * ======================= ===============================================
 *
 * These are wrapper functions to taking steps with several of REBOUND's integrators in order to build custom splitting schemes.
 *
 * **Effect Parameters**
 *
 * None
 * 
 * **Particle Parameters**
 *
 * None
 *
 */

#include <math.h>
#include "rebound.h"
#include "reboundx.h"

// will do IAS with gravity + any additional_forces

void rebx_ias15_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    const double old_t = sim->t;
    const double t_needed = old_t + dt;
    const double old_dt = sim->dt;
    sim->gravity_ignore_terms = 0;
    reb_integrator_ias15_reset(sim);
    
    sim->dt = 0.0001*dt; // start with a small timestep.
    
    while(sim->t < t_needed && fabs(sim->dt/old_dt)>1e-14 ){
        reb_update_acceleration(sim);
        reb_integrator_ias15_part2(sim);
        if (sim->t+sim->dt > t_needed){
            sim->dt = t_needed-sim->t;
        }
    }
    sim->t = old_t;
    sim->dt = old_dt; // reset in case this is part of a chain of steps
}

void rebx_kepler_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    reb_integrator_whfast_init(sim);
    reb_integrator_whfast_from_inertial(sim);
    reb_whfast_kepler_step(sim, dt);
    reb_whfast_com_step(sim, dt);
    reb_integrator_whfast_to_inertial(sim);
}

void rebx_jump_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    reb_integrator_whfast_init(sim);
    reb_integrator_whfast_from_inertial(sim);
    reb_whfast_jump_step(sim, dt);
    reb_integrator_whfast_to_inertial(sim);
}

void rebx_interaction_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    reb_integrator_whfast_init(sim);
    reb_integrator_whfast_from_inertial(sim);
    reb_update_acceleration(sim);
    reb_whfast_interaction_step(sim, dt);
    reb_integrator_whfast_to_inertial(sim);
}

void rebx_drift_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
	const int N = sim->N;
	struct reb_particle* restrict const particles = sim->particles;
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].x  += dt * particles[i].vx;
		particles[i].y  += dt * particles[i].vy;
		particles[i].z  += dt * particles[i].vz;
	}
}

void rebx_kick_step(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    reb_update_acceleration(sim);
	const int N = sim->N;
	struct reb_particle* restrict const particles = sim->particles;
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].vx += dt * particles[i].ax;
		particles[i].vy += dt * particles[i].ay;
		particles[i].vz += dt * particles[i].az;
	}
}
