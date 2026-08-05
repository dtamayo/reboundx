/**
 * @file    radiation_forces.c
 * @brief   Add radiation forces
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
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Radiation Forces$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 H. Rein, D. Tamayo
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_.
 * Based on                `Burns et al. 1979 <http://labs.adsabs.harvard.edu/adsabs/abs/1979Icar...40....1B/>`_.
 * C Example               :ref:`c_example_rad_forces_debris_disk`, :ref:`c_example_rad_forces_circumplanetary`.
 * Python Example          `Radiation_Forces_Debris_Disk.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Radiation_Forces_Debris_Disk.ipynb>`_,
 *                         `Radiation_Forces_Circumplanetary_Dust.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Radiation_Forces_Circumplanetary_Dust.ipynb>`_.
 * ======================= ===============================================
 * 
 * This applies radiation forces to particles in the simulation.  
 * It incorporates both radiation pressure and Poynting-Robertson drag.
 * Only particles whose `beta` parameter is set will feel the radiation.  
 * 
 * **Effect Parameters**
 * 
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * c (double)                   Yes         Speed of light in the units used for the simulation.
 * shadow_model (int)           No          Flag identifying the model for shadows (0 or NULL = no_shadow, 1 = cylindrical, 2 = conical)
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * If no particles have radiation_source set, effect will assume the particle at index 0 in the particles array is the source.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * radiation_source (int)       No          Flag identifying the particle as the source of radiation.
 * beta (float)                 Yes         Ratio of radiation pressure force to gravitational force. Particles without beta set feel no radiation forces.
 * shadow_creator (int)         No          Flag identifying the particle as the source of shadow.
 * R_eq (double)                No          Radius of the body (required for source and occulting bodies if shadowing is used).
 * ============================ =========== ==================================================================
 * 
 * The following changes in the radiation forces aim to establish a "shadow function" in order to achieve a realistic physical approach to the problem. 
 * For achieving this goal, two shadow models are implemented in order for the user to decide which fits best in their simulation, while 
 * preserving backward compatibility with the previous version.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "reboundx.h"


// Helper function to safely degrade to cylindrical model if the source lacks a radius
static int get_active_shadow_flag(struct rebx_extras* const rebx, struct reb_simulation* const sim, int base_flag, struct reb_particle* p, int* warned_degrade) {
    if (base_flag == 2) {
        double* source_req = rebx_get_param(rebx, p->ap, "R_eq");
        if (source_req == NULL || *source_req <= 0.0) {
            if (!(*warned_degrade)) {
                reb_simulation_warning(sim, "Conical shadow requested but radiation_source lacks valid R_eq. Degrading to cylindrical model.");
                *warned_degrade = 1;
            }
            return 1;
        }
    }
    return base_flag;
}




static void rebx_calculate_radiation_forces(struct rebx_extras* const rebx, struct reb_simulation* const sim, 
                                            const double c, const int source_index, struct reb_particle* const particles, 
                                            const int N, const int shadow_flag, const int num_occulters, 
                                            const int* occ_index, const double* occ_radii){

    const struct reb_particle source = particles[source_index];
    const double mu = sim->G*source.m;

    double* source_radius = rebx_get_param(rebx, source.ap, "R_eq");
    const double R_source = (source_radius != NULL) ? *source_radius : 0.0;

    for (int i=0; i<N; i++){
        
        if(i == source_index) continue;
        
        const double* beta = rebx_get_param(rebx, particles[i].ap, "beta");
        if(beta == NULL) continue; // only particles with beta set feel radiation forces

        const double dx = particles[i].x - source.x;
        const double dy = particles[i].y - source.y;
        const double dz = particles[i].z - source.z;
        const double dr = sqrt(dx*dx + dy*dy + dz*dz); // distance to star

        double shadow_factor = 1.0;

        for (int k = 0; k < num_occulters; k++){

            int j = occ_index[k];

            if (j == i || j == source_index ) {
                continue;
            }

            const double R_occ = occ_radii[k];

            const struct reb_particle occulter = particles[j];
            const double occ_x = occulter.x - source.x; 
            const double occ_y = occulter.y - source.y;
            const double occ_z = occulter.z - source.z;
            const double dr_occ = sqrt(occ_x*occ_x + occ_y*occ_y + occ_z*occ_z); 
            
            if (dr_occ <= 0) {
                continue;
            }

            const double ux = occ_x / dr_occ;
            const double uy = occ_y / dr_occ;
            const double uz = occ_z / dr_occ;

            if (shadow_flag == 1){ // Use cylindrical shadow model

                const double s = dx*ux + dy*uy + dz*uz; // Projection of the vector from source to particle onto the line of sight

                if (s <= dr_occ) {
                    continue; 
                }

                // Perperdicular distance from the line of sight to the particle
                const double perp_x = dx - s*ux;
                const double perp_y = dy - s*uy;
                const double perp_z = dz - s*uz;
                const double perp_dist = sqrt(perp_x*perp_x + perp_y*perp_y + perp_z*perp_z);

                if (perp_dist < R_occ) {
                    shadow_factor = 0.0;
                    break;
                }
                
            }

            else if (shadow_flag == 2){ // Use conical shadow model

                const double r_particle_x = particles[i].x - occulter.x;
                const double r_particle_y = particles[i].y - occulter.y;
                const double r_particle_z = particles[i].z - occulter.z;
                const double r_particle_mod = sqrt(r_particle_x*r_particle_x + r_particle_y*r_particle_y + r_particle_z*r_particle_z);

                const double dot = (occ_x*r_particle_x + occ_y*r_particle_y + occ_z*r_particle_z);
                const double s0 = dot/dr_occ; // projection of the vector from occulter to particle onto the line from occulter to source

                if (s0 < 0.0){
                    continue; 
                }

                const double dot_prod = r_particle_x*dx + r_particle_y*dy + r_particle_z*dz;
                const double cross_x = dy*r_particle_z - dz*r_particle_y;
                const double cross_y = dz*r_particle_x - dx*r_particle_z;
                const double cross_z = dx*r_particle_y - dy*r_particle_x;
                const double cross_mag = sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);

                const double app_r_source = asin(fmin(R_source / dr, 1.0));           
                const double app_r_occ = asin(fmin(R_occ / r_particle_mod, 1.0));        

                const double d_ang = atan2(cross_mag, dot_prod);       // Angular distance between both centers

                // Gnomonic projection for flat-circle geometry 
                // Exact at both boundaries; linear interpolation of the offset in between.
                const double p_src = app_r_source;
                const double p_occ = tan(app_r_occ);
                const double p_dang = d_ang - app_r_occ + p_occ;

                if (p_dang >= p_src + p_occ) {
                    continue; 
                }

                else if (p_dang <= fabs(p_src - p_occ)) {
                    if (p_src <= p_occ) {
                        // Total Eclipse (Umbra)
                        shadow_factor = 0.0; 
                        break;
                    } else {
                        // Annular Eclipse
                        const double f_annular = 1.0 - (p_occ * p_occ) / (p_src * p_src);
                        if (f_annular < shadow_factor) shadow_factor = f_annular;
                    }
                } 

                else {
                    // Partial Eclipse (Gloom) - Two circles intersection 
                    const double d1 = (p_src*p_src - p_occ*p_occ + p_dang*p_dang) / (2.0 * p_dang);
                    const double d2 = p_dang - d1;
                    
                    double arg1 = d1 / p_src;
                    if (arg1 > 1.0) arg1 = 1.0; else if (arg1 < -1.0) arg1 = -1.0;
                    
                    double arg2 = d2 / p_occ;
                    if (arg2 > 1.0) arg2 = 1.0; else if (arg2 < -1.0) arg2 = -1.0;

                    // Overlapping areas
                    const double area = p_src * p_src * acos(arg1) - d1 * sqrt(fabs(p_src*p_src - d1*d1))
                                    + p_occ * p_occ * acos(arg2) - d2 * sqrt(fabs(p_occ*p_occ - d2*d2));
                    
                    const double f_partial = 1.0 - area / (M_PI * p_src * p_src);
                    
                    if (f_partial < shadow_factor) {
                        shadow_factor = f_partial;
                    }
                } 
            }            
        }
        
        const double dvx = particles[i].vx - source.vx;
        const double dvy = particles[i].vy - source.vy;
        const double dvz = particles[i].vz - source.vz;
        const double rdot = (dx*dvx + dy*dvy + dz*dvz)/dr; // radial velocity
        const double a_rad = shadow_factor*(*beta)*mu/(dr*dr);

        // Equation (5) of Burns, Lamy & Soter (1979)

        particles[i].ax += a_rad*((1.-rdot/c)*dx/dr - dvx/c);
        particles[i].ay += a_rad*((1.-rdot/c)*dy/dr - dvy/c);
        particles[i].az += a_rad*((1.-rdot/c)*dz/dr - dvz/c);
	}
}

// Hard limit for stack-allocated occulters buffer to avoid VLA
#define MAX_STATIC_OCC 16

void rebx_radiation_forces(struct reb_simulation* const sim, struct rebx_force* const radiation_forces, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    double* c = rebx_get_param(rebx, radiation_forces->ap, "c");
    if (c == NULL){
        reb_simulation_error(sim, "Need to set speed of light in radiation_forces effect.  See examples in documentation.\n");
        return;
    }
    
    int source_found=0;

    int* shadow_flag_ptr = rebx_get_param(rebx, radiation_forces->ap, "shadow_model");
    int shadow_flag = 0;

    // Throttling warnings per condition to avoid stderr spam on every IAS15 substep
    static int warned_model = 0;
    static int warned_occ_req = 0;
    static int warned_degrade = 0;

    // Reset flags at the start of a new simulation to prevent cross-simulation contamination
    if (sim->steps_done == 0 || sim->t == 0.0) {
        warned_model = 0;
        warned_occ_req = 0;
        warned_degrade = 0;
    }

    if (shadow_flag_ptr != NULL) {
        if (*shadow_flag_ptr < 0 || *shadow_flag_ptr > 2){
            if (!warned_model) {
                reb_simulation_warning(sim, "REBOUNDx: shadow_model must be 0, 1 or 2. Shadowing disabled.\n");
                warned_model = 1;
            }        
        } else{
            shadow_flag = *shadow_flag_ptr;
        }            
    }

    int occ_index_static[MAX_STATIC_OCC];
    double occ_radii_static[MAX_STATIC_OCC];
    int* occ_index = occ_index_static;
    double* occ_radii = occ_radii_static;
    
    int num_occulters = 0;
    int allocated_dynamically = 0;

    if (shadow_flag > 0) {        
        for (int i = 0; i < N; i++) {
            int* shadow_creator = rebx_get_param(rebx, particles[i].ap, "shadow_creator");
            if (shadow_creator != NULL && *shadow_creator != 0) {
                double* req = rebx_get_param(rebx, particles[i].ap, "R_eq");
                if (req != NULL) {
                    // Fallback to malloc if we exceed the static stack buffer (very rare)
                    if (num_occulters == MAX_STATIC_OCC) {
                        occ_index = malloc(N * sizeof(int));
                        occ_radii = malloc(N * sizeof(double));
                        if (occ_index == NULL || occ_radii == NULL) {
                            reb_simulation_error(sim, "Memory allocation failed for occulters array. Aborting forces evaluation.");
                            if (occ_index != NULL && occ_index != occ_index_static) free(occ_index);
                            if (occ_radii != NULL && occ_radii != occ_radii_static) free(occ_radii);
                            return;
                        }
                        // Copy existing static data to the newly allocated dynamic array
                        for (int k = 0; k < MAX_STATIC_OCC; k++) {
                            occ_index[k] = occ_index_static[k];
                            occ_radii[k] = occ_radii_static[k];
                        }
                        allocated_dynamically = 1;
                    }
                    occ_index[num_occulters] = i;
                    occ_radii[num_occulters] = *req;
                    num_occulters++;
                } else if (!warned_occ_req) {
                    reb_simulation_warning(sim, "Particle has shadow_creator set but no R_eq defined. It will be ignored as an occulter.");
                    warned_occ_req = 1;
                }
            }
        }
    }


    // Loop for radiation sources
    for (int i = 0; i < N; i++) {
        if (rebx_get_param(rebx, particles[i].ap, "radiation_source") != NULL) {
            source_found = 1;
            int active_shadow_flag = get_active_shadow_flag(rebx, sim, shadow_flag, &particles[i], &warned_degrade);
            rebx_calculate_radiation_forces(rebx, sim, *c, i, particles, N, active_shadow_flag, num_occulters, occ_index, occ_radii);
        }
    }

    // Preserving backward compatibility: Default to index 0 if no radiation_source is set
    if (source_found == 0 && N > 0) {
        int active_shadow_flag = get_active_shadow_flag(rebx, sim, shadow_flag, &particles[0], &warned_degrade);
        rebx_calculate_radiation_forces(rebx, sim, *c, 0, particles, N, active_shadow_flag, num_occulters, occ_index, occ_radii);
    }

    if (allocated_dynamically) {
        free(occ_index);
        free(occ_radii);
    }
}

double rebx_rad_calc_beta(const double G, const double c, const double source_mass, const double source_luminosity, const double radius, const double density, const double Q_pr){
    return 3.*source_luminosity*Q_pr/(16.*M_PI*G*source_mass*c*density*radius);   
}
double rebx_rad_calc_particle_radius(const double G, const double c, const double source_mass, const double source_luminosity, const double beta, const double density, const double Q_pr){
    return 3.*source_luminosity*Q_pr/(16.*M_PI*G*source_mass*c*density*beta);
}