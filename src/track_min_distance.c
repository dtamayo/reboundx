/**
 * @file    track_min_distance.c
 * @brief   Track minimum distance of secondaries from primary each timestep and log results.
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
#include "rebound.h"
#include "reboundx.h"

void rebx_track_min_distance(struct reb_simulation* const sim, struct rebx_effect* const effect){
    const int _N_real = sim->N - sim->N_var;	
	for(int i=0; i<_N_real; i++){
		struct reb_particle* const p = &sim->particles[i];
        double* min_distance = rebx_get_param_check(p, "min_distance", REBX_TYPE_DOUBLE);
        if (min_distance != NULL){
		    const uint32_t* const target = rebx_get_param_check(p, "min_distance_from", REBX_TYPE_UINT32);
            struct reb_particle* source;
            if (target == NULL){
                source = &sim->particles[0];
            }
            else{
                source = reb_get_particle_by_hash(sim, *target);
            }
            const double dx = p->x-source->x;
            const double dy = p->y-source->y;
            const double dz = p->z-source->z;
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < *min_distance*(*min_distance)){
                *min_distance = sqrt(r2);
            }
        }
	}
}

