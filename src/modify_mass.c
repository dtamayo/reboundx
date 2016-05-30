/**
 * @file    modify_mass.c
 * @brief   Add exponential mass loss/growth between timesteps in the simulation.
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
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
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "modify_mass.h"
#include "rebound.h"
#include "reboundx.h"

#define TWOPI 6.2831853071795862

void rebx_add_modify_mass(struct rebx_extras* rebx){
    rebx_add_post_timestep_modification(rebx, NULL, "modify_mass", rebx_modify_mass);
}

void rebx_modify_mass(struct reb_simulation* const sim, struct rebx_effect* const effect){
    const int _N_real = sim->N - sim->N_var;	
	const double dt = sim->dt_last_done;
	for(int i=0; i<_N_real; i++){
		struct reb_particle* const p = &sim->particles[i];
        const double tau_mass = rebx_get_param_double(p, "tau_mass");
	    if(isnan(tau_mass)) continue; // only particles with tau_mass set experience mass loss/growth

		p->m += p->m*dt/tau_mass;
	}
    reb_move_to_com(sim);
}

