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
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Miscellaneous Utilities$     // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script). 
 * 
 * ======================= ===============================================
 * Authors                 D. Tamayo
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
 * Based on                None
 * C Example               :ref:`c_example_track_min_distance`
 * Python Example          `TrackMinDistance.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TrackMinDistance.ipynb>`_.
 * ======================= ===============================================
 * 
 * For a given particle, this keeps track of that particle's minimum distance from another body in the simulation.  User
 * should add parameters to the particular particle whose distance should be tracked.
 *
 * **Effect Parameters**
 * 
 * *None*
 * 
 * **Particle Parameters**
 * 
 * Only particles with their ``min_distance`` parameter set initially will track their minimum distance. The effect will
 * update this parameter when the particle gets closer than the value of ``min_distance``, so the user has to set it
 * initially.  By default, distance is measured from sim->particles[0], but you can specify a different particle by setting
 * the ``min_distance_from`` parameter to the hash of the target particle.
 * 
 * ================================ =========== =======================================================
 * Name (C type)                    Required    Description
 * ================================ =========== =======================================================
 * min_distance (double)            Yes         Particle's mininimum distance.
 * min_distance_from (uint32)       No          Hash for particle from which to measure distance
 * min_distance_orbit (reb_orbit)   No          Parameter to store orbital elements at moment corresponding to min_distance (heliocentric)
 * ================================ =========== =======================================================
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_track_min_distance(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    const struct rebx_extras* const rebx = sim->extras;
    const int N = sim->N - sim->N_var;
    for(int i=0; i<N; i++){
        struct reb_particle* const p = &sim->particles[i];
        double* min_distance = rebx_get_param(rebx, p->ap, "min_distance");
        if (min_distance != NULL){
            const uint32_t* const target = rebx_get_param(rebx, p->ap, "min_distance_from");
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
                struct reb_orbit* const orbit = rebx_get_param(rebx, p->ap, "min_distance_orbit");
                if (orbit != NULL){
                    *orbit = reb_tools_particle_to_orbit(sim->G, *p, *source);
                }
            }
        }
    }
}

