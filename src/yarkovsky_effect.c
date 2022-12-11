/**
 * @file    yarkovsky_effect.c
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
 * $Radiation Forces$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 Noah Ferich, D. Tamayo
 * Implementation Paper    Ferich et al., in prep.
 * Based on                `Veras et al., 2015 <https://ui.adsabs.harvard.edu/abs/2015MNRAS.451.2814V/abstract>`_, `Veras et al., 2019 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.485..708V/abstract>`_.
 * C Example               :ref:`c_example_yarkovsky_effect`.
 * Python Example          `YarkovskyEffect.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/YarkovskyEffect.ipynb>`_.
 * ======================= ===============================================
 *
 * Adds the accelerations and orbital perturbations created by the Yarkovsky effect onto one or more bodies in the simulation. There are two distinct versions of this effect that can be used: the 'full version' and the 'simple version'. The full version uses the full equations found in Veras et al. (2015) to accurately calculate the Yarkovsky effect on a particle. However, this version slows down simulations and requies a large amount of parameters. For these reasons, the simple version of the effect (based on Veras et al. (2019)) is available. While the magnitude of the acceleration created by the effect will be the same, this version places constant values in a crucial rotation matrix to simplify the push from the Yarkovsky effect on a body. This version is faster and requires less parameters and can be used to get an upper bound on how much the Yarkovsky effect can push an object's orbit inwards or outwards. The lists below describes which parameters are needed for one or both versions of this effect. For more information, please visit the papers and examples linked above.
 *
 * **Effect Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * ye_lstar (float)             Yes         Luminosity of sim's star (Required for both versions).
 * ye_c (float)                 Yes         Speed of light (Required for both versions).
 * ye_stef_boltz (float)        No          Stefan-Boltzmann constant (Required for full version).
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * particles[i].r (float)       Yes         Physical radius of a body (Required for both versions).
 * ye_flag (int)                Yes         0 sets full version of effect. 1 uses simple version with outward migration. -1 uses the simple version with inward migration (see examples and paper).
 * ye_body_density (float)      Yes         Density of an object (Required for both versions)
 * ye_rotation_period (float)   No          Rotation period of a spinning object (Required for full version)
 * ye_albedo (float)            Yes         Albedo of an object (Reuired for both versions)
 * ye_emissivity (float)        No          Emissivity of an object (Required for full version)
 * ye_thermal_inertia (float)   No          Thermal inertia of an object (Required for full version)
 * ye_k (float)                 No          A constant that gets a value between 0 and 1/4 based on the object's rotation - see Veras et al. (2015) for more information on it (Required for full version)
 * ye_spin_axis_x (float)       No          The x value for the spin axis vector of an object (Required for full version)
 * ye_spin_axis_y (float)       No          The y value for the spin axis vector of an object (Required for full version)
 * ye_spin_axis_z (float)       No          The z value for the spin axis vector of an object (Required for full version)
 * ============================ =========== ==================================================================
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "reboundx.h"

static void rebx_calculate_yarkovsky_effect(struct reb_simulation* sim, struct reb_particle* target, struct reb_particle* star, double G, double *density, double *lstar, double *rotation_period, double *Gamma, double *albedo, double *emissivity, double *k, double *c, double *stef_boltz, int *yark_flag, double *sx, double *sy, double *sz){
    
    int i; //variables needed for future iteration loops
    int j;
    double unit_matrix[3][3] = {{1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{0.0, 0.0, 1.0}};

    double radius = target->r;
    
    double q_yar = 1.0-(*albedo);
    
    double dx = target->x - star->x;
    double dy = target->y - star->y;
    double dz = target->z - star->z;
    
    double dvx = target->vx - star->vx;
    double dvy = target->vy - star->vy;
    double dvz = target->vz - star->vz;
    
    double distance = sqrt((dx*dx)+(dy*dy)+(dz*dz)); //distance of asteroid from the star
    
    double rdotv = ((dx*dvx)+(dy*dvy)+(dz*dvz))/((*c)*distance); //dot product of position and velocity vectors- the term in the denominator is needed when calculating the i-vector
    
    double i_vector[3][1];
    
    i_vector[0][0] = ((1-rdotv)*(dx/distance))-(dvx/(*c));
    i_vector[1][0] = ((1-rdotv)*(dy/distance))-(dvy/(*c));
    i_vector[2][0] = ((1-rdotv)*(dz/distance))-(dvz/(*c));
    
    double yarkovsky_magnitude; //magnitude of force created by the effect

    double yark_matrix[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
    
    
    if (*yark_flag == 1) {
        
        yark_matrix[1][0] = 1.0; // maximizes the effect pushing outwards
        
        yarkovsky_magnitude = (3*q_yar*(*lstar))/(64*M_PI*radius*(*density)*(*c)*distance*distance);
    }
    
    if (*yark_flag == -1) {
        
        yark_matrix[0][1] = 1.0; //maximizes the effect pushing inwards
        
        yarkovsky_magnitude = (3*q_yar*(*lstar))/(64*M_PI*radius*(*density)*(*c)*distance*distance);
    }
    
    //will run through full equations to create the yark_matrix
    if (*yark_flag == 0) {
        
        //makes sure all necessary parameters have been entered
        if (stef_boltz == NULL || rotation_period == NULL || Gamma == NULL || albedo == NULL || emissivity == NULL || k == NULL || sx == NULL || sy == NULL || sz == NULL) {
            reb_error(sim, "REBOUNDx Error: One or more parameters missing for this version of the Yarkovsky effect in Rebx. Please make sure you've given values to all variables for this version before running simulations. See documentation and YarkovskyEffect.ipynb. If you'd rather use the simplified version of this effect (requires fewer parameters), then please set 'yark_flag' to -1 or 1.\n\n");
            return;
        }
        
        struct reb_orbit o = reb_tools_particle_to_orbit(G, *target, *star);
        
        yarkovsky_magnitude = (3*(*k)*q_yar*(*lstar))/(16*M_PI*radius*(*density)*(*c)*distance*distance);

        double Smag = sqrt(((*sx)*(*sx))+ (*sy)*(*sy) + (*sz)*(*sz));
    
        double hx = (dy*dvz)-(dz*dvy);
        double hy = (dz*dvx)-(dx*dvz);
        double hz = (dx*dvy)-(dy*dvx);
        double Hmag = sqrt((hx*hx)+ (hy*hy) + (hz*hz));
    
        double inv_smag = 1.0/Smag;
        double inv_mag_sqrd = 1.0/(Smag*Smag);
        double inv_hmag = 1.0/Hmag;
        double inv_hmag_sqrd = 1.0/(Hmag*Hmag);
    
        double R1s[3][3] = {{0.0, -(*sz)*inv_smag, (*sy)*inv_smag},{(*sz)*inv_smag, 0.0, -(*sx)*inv_smag},{-(*sy)*inv_smag, (*sx)*inv_smag, 0.0}};
        
        double R2s[3][3] = {{(*sx)*(*sx)*inv_mag_sqrd, (*sx)*(*sy)*inv_mag_sqrd, (*sx)*(*sz)*inv_mag_sqrd},{(*sx)*(*sy)*inv_mag_sqrd, (*sy)*(*sy)*inv_mag_sqrd, (*sy)*(*sz)*inv_mag_sqrd},{(*sx)*(*sz)*inv_mag_sqrd, (*sy)*(*sz)*inv_mag_sqrd, (*sz)*(*sz)*inv_mag_sqrd}};
        
        double R1h[3][3] = {{0.0, -hz*inv_hmag, hy*inv_hmag},{hz*inv_hmag, 0.0, -hx*inv_hmag},{-hy*inv_hmag, hx*inv_hmag, 0.0}};
        
        double R2h[3][3] = {{hx*hx*inv_hmag_sqrd, hx*hy*inv_hmag_sqrd, hx*hz*inv_hmag_sqrd},{hx*hy*inv_hmag_sqrd, hy*hy*inv_hmag_sqrd, hy*hz*inv_hmag_sqrd},{hx*hz*inv_hmag_sqrd, hy*hz*inv_hmag_sqrd, hz*hz*inv_hmag_sqrd}};

        double tanPhi = 1.0/(1.0+(.5*pow(((*stef_boltz)*(*emissivity))/(M_PI*M_PI*M_PI*M_PI*M_PI), .25))*sqrt((*rotation_period)/((*Gamma)*(*Gamma)))*pow((*lstar*q_yar)/(distance*distance), .75));
    
        double tanEpsilon = 1.0/(1.0+(.5*pow((*stef_boltz*(*emissivity))/(M_PI*M_PI*M_PI*M_PI*M_PI), .25))*sqrt((o.P)/((*Gamma)*(*Gamma)))*pow((*lstar*q_yar)/(distance*distance), .75));
    
        double Phi = atan(tanPhi);
        double Epsilon = atan(tanEpsilon);
        
        double cos_phi = cos(Phi);
        double sin_phi = sin(Phi);
        double cos_epsilon = cos(Epsilon);
        double sin_epsilon = sin(Epsilon);
    
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
                yark_matrix[i][j] = (Rys[i][0]*Ryh[0][j]) + (Rys[i][1]*Ryh[1][j]) + (Rys[i][2]*Ryh[2][j]);
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
        
        double* density = rebx_get_param(rebx, target->ap, "ye_body_density");
        double* lstar = rebx_get_param(rebx, force->ap, "ye_lstar");
        double* rotation_period = rebx_get_param(rebx, target->ap, "ye_rotation_period");
        double* Gamma = rebx_get_param(rebx, target->ap, "ye_thermal_inertia");
        double* albedo = rebx_get_param(rebx, target->ap, "ye_albedo");
        double* emissivity = rebx_get_param(rebx, target->ap, "ye_emissivity");
        double* k = rebx_get_param(rebx, target->ap, "ye_k");
        double* c = rebx_get_param(rebx, force->ap, "ye_c");
        double* stef_boltz = rebx_get_param(rebx, force->ap, "ye_stef_boltz");
        int* yark_flag = rebx_get_param(rebx, target->ap, "ye_flag");
        double* sx = rebx_get_param(rebx, target->ap, "ye_spin_axis_x");
        double* sy = rebx_get_param(rebx, target->ap, "ye_spin_axis_y");
        double* sz = rebx_get_param(rebx, target->ap, "ye_spin_axis_z");
        
        //if these necessary conditions are met the Yarkovsky effect will be calculated for a particle in the sim
        if (density != NULL && target->r != 0 && albedo != NULL && lstar != NULL && c != NULL && yark_flag != NULL){
            rebx_calculate_yarkovsky_effect(sim, target, star, G, density, lstar, rotation_period, Gamma, albedo, emissivity, k, c, stef_boltz, yark_flag, sx, sy, sz);
        }
    }
}
