/**
 *
 *
 * @file    tides_constant_time_lag.c
 * @brief   Add constant time lag tides raised on primary, orbiting bodies, or both
 * @author  Stanley A. Baronett <stanley.a.baronett@gmail.com>, Dan Tamayo <tamayo.daniel@gmail.com>
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
 * $Tides$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 Stanley A. Baronett, D. Tamayo, Noah Ferich
 * Implementation Paper    Baronett et al., in prep.
 * Based on                `Hut 1981 <https://ui.adsabs.harvard.edu/#abs/1981A&A....99..126H/abstract>`_, `Bolmont et al., 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...583A.116B/abstract>`_.
 * C Example               :ref:`c_example_tides_constant_time_lag`.
 * Python Example          `TidesConstantTimeLag.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesConstantTimeLag.ipynb>`_.
 * ======================= ===============================================
 *
 * This adds constant time lag tidal interactions between orbiting bodies in the simulation and the primary, both from tides raised on the primary and on the other bodies.
 * In all cases, we need to set masses for all the particles that will feel these tidal forces. After that, we can choose to include tides raised on the primary, on the "planets", or both, by setting the respective bodies' physical radius particles[i].r, k1 (apsidal motion constant, half the tidal Love number), constant time lag tau, and rotation rate Omega. See Hut (1981) and Bolmont et al. 2015 above.
 *
 * If tau is not set, it will default to zero and yield the conservative piece of the tidal potential.
 * 
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * particles[i].r (float)       Yes         Physical radius (required for contribution from tides raised on the body).
 * tctl_k1 (float)                   Yes         Apsidal motion constant (half the tidal Love number k2).
 * tctl_tau (float)                  No          Constant time lag. If not set will default to 0 and give conservative tidal potential
 * Omega (float)                No          Rotation rate. If not set will default to 0
 * ============================ =========== ==================================================================
 *
 *
 * CHANGE THIS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "reboundx.h"

static void rebx_calculate_max_yarkovsky_effect(struct reb_particle* target, double G, double *density, double *lstar, double *c, int *direction_flag){
    
    double yark_matrix[3][3] = {{1, 0, 0},{0, 1, 0},{0, 0, 1}};
    
    if (*direction_flag == 1) {
        // outward drift matrix is {{1, 0, 0},{.25, 1, 0},{0, 0, 1}};
        
        yark_matrix[1][0] = .25;
    }
    
    if (*direction_flag == -1){
        // inward drift matrix is {{1, .25, 0},{0, 1, 0},{0, 0, 1}};
        
        yark_matrix[0][1] = .25;
    }
    
    int i; //variable needed for future iteration loops
    
    double radius = target->r;
    
    double distance = sqrt((target->x*target->x)+(target->y*target->y)+(target->z*target->z)); //distance of asteroid from the star
    
    double rdotv = ((target->x*target->vx)+(target->y*target->vy)+(target->z*target->vz))/((*c)*distance); //dot product of position and velocity vectors- the term in the denominator is needed when calculating the i vector
    
    double i_vector[3][1];
    
    i_vector[0][0] = ((1-rdotv)*(target->x/distance))-(target->vx/(*c));
    i_vector[1][0] = ((1-rdotv)*(target->y/distance))-(target->vy/(*c));
    i_vector[2][0] = ((1-rdotv)*(target->z/distance))-(target->vz/(*c));
    
    double yarkovsky_magnitude = (*lstar)/(4*((4*M_PI*radius*(*density))/3)*(*c)*distance*distance);
    
    double direction_matrix[3][1];
    
    //loops calcuates a vector which gives the direction of the acceleration created by the Yarkovsky effect
    for (i=0; i<3; i++){
        direction_matrix[i][0] = (yark_matrix[i][0]*i_vector[0][0]) + (yark_matrix[i][1]*i_vector[1][0]) + (yark_matrix[i][2]*i_vector[2][0]);
    }
    
    double yarkovsky_acceleration[3][1];
    
    for (i=0; i<3; i++){
     
        yarkovsky_acceleration[i][0] = (direction_matrix[i][0]*yarkovsky_magnitude);
        
    }
    
        //adds Yarkovsky aceleration to the asteroid's acceleration in the sim
        target->ax += yarkovsky_acceleration[0][0];
        target->ay += yarkovsky_acceleration[1][0];
        target->az += yarkovsky_acceleration[2][0];

    }

void rebx_max_yarkovsky(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
        
    struct rebx_extras* const rebx = sim->extras;
    double G = sim->G;
    
    for (int i=1; i<N; i++){
        struct reb_particle* target = &particles[i];
        double* density = rebx_get_param(rebx, target->ap, "my_body_density");
        double* lstar = rebx_get_param(rebx, force->ap, "my_lstar");
        double* c = rebx_get_param(rebx, force->ap, "my_c");
        int* direction_flag = rebx_get_param(rebx, target->ap, "direction_flag");
        if (density != NULL && target->r != 0 && lstar != NULL && direction_flag != NULL && c != NULL){
            rebx_calculate_max_yarkovsky_effect(target, G, density, lstar, c, direction_flag);
        }
        
    }
}
