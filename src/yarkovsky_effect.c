//
//  rebx_yarkovsky.c
//  
//
//  Created by Noah Ferich on 5/17/21.
//

//#include "rebx_yarkovsky.h"

/**
 *
 *
 * @file    yarkovsky.c
 * @brief   Adds the perturbations from the Yarkovsky effect to one or more of the orbiting bodies
 * @author  Noah Ferich <nofe4108@colorado.edu>, Dan Tamayo <tamayo.daniel@gmail.com>
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
 * Based on                `Veras et al., 2015 <https://ui.adsabs.harvard.edu/abs/2015MNRAS.451.2814V/abstract>`_, `Veras et al., 2019 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.485..708V/abstract>`_.
 * C Example               :ref:`c_example_tides_constant_time_lag`.
 * Python Example          `TidesConstantTimeLag.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesConstantTimeLag.ipynb>`_.
 * ======================= ===============================================
 *
 * This adds the accelerations and orbital changes created by the Yarkovsky effect onto one or more bodies in the simulation.
 * There are two distinct versions of this effect that can be used. One version uses the full equations found in Veras et. al (2015) to accurately calculate the Yarkovsky effect on a particle. However, this version slows down simulations and requies a large amount of parameters. For these reasons, a second version of the effect is available and is described in Veras et. al (2019). While the magnitude of the acceleration created by the effect will be the same, this version artifically places values in a crucial rotation matrix to maximize the push from the Yarkovsky effect on a body. This version is faster and requires less parameters and can be used to get an upper bound on how much the Yarkovsky effect can push an object's orbit inwards or outwards. The list below describes which parameters are needed for one of both versions of this effect. For more information, please visit the papers linked above.
 *
 * As is standard with all Rebx effects, the parameters for an effect must be entered using the same units as the simulation.
 *
 * **Effect Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * lstar (float)                   Yes              Luminosity of sim's star (Required for both versions).
 * yark_c (float)                Yes             Speed of light (Required for both versions).
 * ============================ =========== ==================================================================
 *
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * particles[i].r (float)       Yes         Physical radius (Required for both versions).
 *
 * ============================ =========== ==================================================================
 *
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "reboundx.h"

static void rebx_calculate_yarkovsky_effect(struct reb_particle* target, struct reb_particle* star, double G, double *density, double *lstar, double *rotation_period, double *C, double *K, double *albedo, double *emissivity, double *k, double *c, double *stef_boltz, int *yark_flag){
    
    int i; //variable needed for future iteration loops
    int j;
    double unit_matrix[3][3] = {{1.0, 1.0, 1.0},{1.0, 1.0, 1.0},{1.0, 1.0, 1.0}};

    double radius = target->r;
    
     double distance = sqrt((target->x*target->x)+(target->y*target->y)+(target->z*target->z)); //distance of asteroid from the star
    
    double rdotv = ((target->x*target->vx)+(target->y*target->vy)+(target->z*target->vz))/((*c)*distance); //dot product of position and velocity vectors- the term in the denominator is needed when calculating the i-vector
    
    double i_vector[3][1];
    
    i_vector[0][0] = ((1-rdotv)*(target->x/distance))-(target->vx/(*c));
    i_vector[1][0] = ((1-rdotv)*(target->y/distance))-(target->vy/(*c));
    i_vector[2][0] = ((1-rdotv)*(target->z/distance))-(target->vz/(*c));
    
    double yarkovsky_magnitude = (*lstar)/(4*((4*M_PI*radius*(*density))/3)*(*c)*distance*distance); //magnitude of force created by the effect

    double yark_matrix[3][3] = {{1, 0, 0},{0, 1, 0},{0, 0, 1}};
    
    
    if (*yark_flag == 1) {
        
        yark_matrix[1][0] = .25; // maximizes the effect pushing outwards
    }
    
    if (*yark_flag == -1) {
        
        yark_matrix[0][1] = .25; //maximizes the effect pushing inwards
    }
    
    //will run through full equations to create the yark_matrix
    if (*yark_flag == 0) {
        
        if (stef_boltz == NULL || rotation_period == NULL || C == NULL || K == NULL || albedo == NULL || emissivity == NULL || k == NULL) {
            
            printf("ERROR: One or more parameters missing for this version of the Yarkovsky effect in Rebx. Please make sure you've given values to all variables for this version before running simulations. If you'd rather use the simplified version of this effect (requires fewer parameters), then please set 'yark_flag' to -1 or 1.\n\n");
            
            printf("MISSING PARAMETERS: \n");
            
            if (stef_boltz == NULL){
                printf("stef_boltz\n");
            }
            if (rotation_period == NULL){
                printf("rotation_period\n");
            }
            if (C == NULL){
                printf("specific_heat_capacity\n");
            }
            if (K == NULL){
                printf("thermal_conductivity\n");
            }
            if (albedo == NULL){
                printf("albedo\n");
            }
            if (emissivity == NULL){
                printf("emissivity\n");
            }
            if (k == NULL){
                printf("k\n");
            }
            exit(0);
        }
        
        struct reb_orbit o = reb_tools_particle_to_orbit(G, *target, *star);
        
        
        double q_yar = 1.0-(*albedo);
    
        double sx = 0.0872;
        double sy = 0.0;
        double sz = -0.9962;
        double Smag = sqrt((sx*sx)+ (sy*sy) + (sz*sz));
    
        double hx = (target->y*target->vz)-(target->z*target->vy);
        double hy = (target->z*target->vx)-(target->x*target->vz);
        double hz = (target->x*target->vy)-(target->y*target->vx);
        double Hmag = sqrt((hx*hx)+ (hy*hy) + (hz*hz));
    
        double inv_smag = 1.0/Smag;
        double inv_mag_sqrd = 1.0/(Smag*Smag);
        double inv_hmag = 1.0/Hmag;
        double inv_hmag_sqrd = 1.0/(Hmag*Hmag);
    
        double R1s[3][3] = {{0.0, -sz*inv_smag, sy*inv_smag},{sz*inv_smag, 0.0, -sx*inv_smag},{-sy*inv_smag, sx*inv_smag, 0.0}};
        double R2s[3][3] = {{sx*sx*inv_mag_sqrd, sx*sy*inv_mag_sqrd, sx*sz*inv_mag_sqrd},{sx*sy*inv_mag_sqrd, sy*sy*inv_mag_sqrd, sy*sz*inv_mag_sqrd},{sx*sz*inv_mag_sqrd, sy*sz*inv_mag_sqrd, sz*sz*inv_mag_sqrd}};
        double R1h[3][3] = {{0.0, -hz*inv_hmag, hy*inv_hmag},{hz*inv_hmag, 0.0, -hx*inv_hmag},{-hy*inv_hmag, hx*inv_hmag, 0.0}};
        double R2h[3][3] = {{hx*hx*inv_hmag_sqrd, hx*hy*inv_hmag_sqrd, hx*hz*inv_hmag_sqrd},{hx*hy*inv_hmag_sqrd, hy*hy*inv_hmag_sqrd, hy*hz*inv_hmag_sqrd},{hx*hz*inv_hmag_sqrd, hy*hz*inv_hmag_sqrd, hz*hz*inv_hmag_sqrd}};
    
        double tanPhi = 1.0/(1.0+(.5*pow(((*stef_boltz)*(*emissivity))/(M_PI*M_PI*M_PI*M_PI*M_PI), .25))*sqrt((*rotation_period)/((*C)*(*K)*(*density)))*pow((*lstar*q_yar)/(distance*distance), .75));
    
        double tanEpsilon = 1.0/(1.0+(.5*pow((*stef_boltz*(*emissivity))/(M_PI*M_PI*M_PI*M_PI*M_PI), .25))*sqrt((o.P)/((*C)*(*K)*(*density)))*pow((*lstar*q_yar)/(distance*distance), .75));
    
        double Phi = atan(tanPhi);
        double Epsilon = atan(tanEpsilon);
    
        double cos_phi = cos(Phi);
        double sin_phi = sin(Phi);
        double sin_epsilon = sin(Epsilon);
        double cos_epsilon = cos(Epsilon);
    
        double Rys[3][3]; //diurnal conntribution for effect
    
        for (i=0; i<3; i++){
            for (j=0; j<3; j++){
                Rys[i][j] = (cos_phi*unit_matrix[i][j]) + (sin_phi*R1s[i][j]) + ((1.0-cos_phi)*R2s[i][j]);
            }
        }
    
        double Ryh[3][3];
    
        for (i=0; i<3; i++){ //seasonal contribution for effect
            for (j=0; j<3; j++){
                Ryh[i][j] = (cos_epsilon*unit_matrix[i][j]) - (sin_epsilon*R1h[i][j]) + ((1-cos_epsilon)*R2h[i][j]);
            }
        }
    
        for (i=0; i<3; i++){
            for(j=0; j<3; j++){
                yark_matrix[i][j] = (Ryh[i][0]*Rys[0][j]) + (Ryh[i][1]*Rys[1][j]) + (Ryh[i][2]*Rys[2][j]);
            }
        }
        
    }
    
    double direction_matrix[3][1];
    
    //calcuates a vector which gives the direction of the acceleration created by the effect
    for (i=0; i<3; i++){
        direction_matrix[i][0] = (yark_matrix[i][0]*i_vector[0][0]) + (yark_matrix[i][1]*i_vector[1][0]) + (yark_matrix[i][2]*i_vector[2][0]);
    }
    
    double yarkovsky_acceleration[3][1];
    
    for (i=0; i<3; i++){
     
        //final result for particle's change in acceleration due to the effect
        yarkovsky_acceleration[i][0] = (yarkovsky_magnitude*direction_matrix[i][0]);
        
    }
    
        //adds Yarkovsky aceleration to the asteroid's acceleration in the sim
        target->ax += yarkovsky_acceleration[0][0];
        target->ay += yarkovsky_acceleration[1][0];
        target->az += yarkovsky_acceleration[2][0];
    
    }

void rebx_yarkovsky_effect(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
        
    struct rebx_extras* const rebx = sim->extras;
    double G = sim->G;
    
    
    for (int i=1; i<N; i++){
        
        struct reb_particle* target = &particles[i];
        struct reb_particle* star = &particles[0];
        
        double* density = rebx_get_param(rebx, target->ap, "body_density");
        double* lstar = rebx_get_param(rebx, force->ap, "lstar");
        double* rotation_period = rebx_get_param(rebx, target->ap, "rotation_period");
        double* C = rebx_get_param(rebx, target->ap, "specific_heat_capacity");
        double* K = rebx_get_param(rebx, target->ap, "thermal_conductivity");
        double* albedo = rebx_get_param(rebx, target->ap, "albedo");
        double* emissivity = rebx_get_param(rebx, target->ap, "emissivity");
        double* k = rebx_get_param(rebx, target->ap, "k");
        double* c = rebx_get_param(rebx, force->ap, "yark_c");
        int* yark_flag = rebx_get_param(rebx, target->ap, "yark_flag");
        double* stef_boltz = rebx_get_param(rebx, force->ap, "stef_boltz");
        
        //if these necessary conditions are met the Yarkovsky effect will be calculated for a particle in the sim
        if (density != NULL && target->r != 0 && lstar != NULL && c != NULL && yark_flag != NULL){
            rebx_calculate_yarkovsky_effect(target, star, G, density, lstar, rotation_period, C, K, albedo, emissivity, k, c, stef_boltz, yark_flag);
        }
        
    }
}
