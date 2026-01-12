/**
 * @file    fragmenting_collisions.c
 * @brief   Adds a REBOUNDx collision module which resolves collisions based on their impact energy
 * @author  Haniyeh Tajer <tajer.1@osu.edu>
 * 
 * @section     LICENSE
 * Copyright (c) 2025 Haniyeh Tajer, Tiger Lu, Dan Tamayo, Hanno Rein
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
 * $Collisions$     // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script). 
 * 
 * ======================= ===============================================
 * Authors                 H. Tajer, H. Rein, T. Lu
 * Implementation Paper    `In prep.`_.
 * Based on                `Leinhardt and Stewart 2012 <https://iopscience.iop.org/article/10.1088/0004-637X/745/1/79>`_, `Chambers 2013 <https://www.sciencedirect.com/science/article/pii/S0019103513000754?via%3Dihub>`_.
 * C Example               :ref:`c_example_fragmenting_collisions_embryo_disk`
 * Python Example          `FragmentingCollisions.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/FragmentingCollisions.ipynb>`_.
 * ======================= ===============================================
 * 
 * This module is based on the collision prescription of Leinhardt & Stewart (2012), from now on called LS2012. 
 * After a collision between two bodies is detected, the outcome will depend on the velocity of the collision, mass of the bodies, and the impact angle.
 * Possible outcomes are:
 * 
 * 1. "merging": particles will merge 
 * 
 * 2. "fragmentation": erosion of the target or accretion of projectile material into target (fragmentation)
 * 
 * 3. "elastic bounce": particles will leave each other intact
 * 
 * For more details, refer to the implementation paper (in prep).
 * 
 * **Effect Parameters**
 * 
 * ======================================== =========== ==============================================================================
 * Field (C type)                           Required    Description
 * ======================================== =========== ==============================================================================
 * fc_min_frag_mass (double)                 Yes         Minimum fragment mass allowed in the simulation. 
 * fc_separation_distance_scale (double)     No          Ratio of distance between the COM of newly added fragment and COM of target and projectile over sum of radii of target and projectile
 * fc_rho1 (double)                          No          Density rho_1 used in computing the spherical radius of combined mass of target and projectile if they had density rho_1. Used in computing the disruption criteria for objects with different bulk densities and mass ratios. For more information, refer to Leinhardt and Stewart 2012.
 * fc_cstar (double)                         No          C* is a dimentionless parameter used as a measure of the dissipation of energy within the target. For more information refer to Leinhardt and Stewart 2012.
 * fc_particle_list_file (string)            No          Name of the output file. Example is: "family_tree.csv". If you wish to produce the output file, you need to set this parameter. 
 * fc_id_max (int)                           No          Maximum particle ID assigned so far. This parameter gets updated everytime a new particle ID is assigned.
 * ======================================== =========== ==============================================================================
 * 
 * **Particle Parameters**
 * 
 * ======================================== =========== ==============================================================================
 * Field (C type)                           Required    Description
 * ======================================== =========== ==============================================================================
 * fc_id (int)                              No          Unique particle ID. Everytime a collision happens, new particle IDs are assigned to all the bodies after the collision, and previous IDs are discarded.
 * ======================================== =========== ==============================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include <stdbool.h>

#define MIN(a, b) ((a) > (b) ? (b) : (a))    // Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    // Returns the maximum of a and b

enum COLLISION_TYPE {
    COLLISION_TYPE_MERGE = 1,
    COLLISION_TYPE_MERGE_B = 2,
    COLLISION_TYPE_MERGE_C = 3,
    COLLISION_TYPE_ACCRETION = 4, // Errosion? Meaning unclear.
    COLLISION_TYPE_SUPERCATASTROPHIC = 5,
    COLLISION_TYPE_BOUNCE_F = 6,
    COLLISION_TYPE_GRAZING_G = 7,
    COLLISION_TYPE_BOUNCE_H = 8,
    COLLISION_TYPE_BOUNCE_I = 9,
    COLLISION_TYPE_HITANDRUN = 10,
};

// Helper function to get radius from mass and density
static double get_radii(double m, double rho){
    return pow((3.*m)/(4.*M_PI*rho),1./3.);
}   

// Function to set a new ID for a new particle
int rebx_fragmenting_collisions_set_new_id(struct reb_simulation* sim, struct rebx_collision_resolve* const collision_resolve, struct reb_particle* p){
    int* fc_id_max = rebx_get_param(sim->extras, collision_resolve->ap, "fc_id_max");
    if (!fc_id_max){ // First call? 
        rebx_set_param_int(sim->extras, &collision_resolve->ap, "fc_id_max", 0);
        fc_id_max = rebx_get_param(sim->extras, collision_resolve->ap, "fc_id_max");
    }
    int new_id = *fc_id_max;
    rebx_set_param_int(sim->extras,  (struct rebx_node**) &(p->ap), "fc_id", new_id);
    (*fc_id_max)++;
    return new_id;
}


static void set_fc_ids(struct reb_simulation* sim, struct rebx_collision_resolve* const collision_resolve){
    int id_count = 0;
    int* fc_id_max = rebx_get_param(sim->extras, collision_resolve->ap, "fc_id_max");
    if (!fc_id_max){ // First call? 
        rebx_set_param_int(sim->extras, &collision_resolve->ap, "fc_id_max", id_count);
        fc_id_max = rebx_get_param(sim->extras, collision_resolve->ap, "fc_id_max");
    }
    else{
        id_count = *(fc_id_max) + 1;
    }
    for(int i=0; i<sim->N; i++){    
        struct reb_particle* p = &(sim->particles[i]);
        int* p_ptr = (int*) rebx_get_param(sim->extras, p->ap, "fc_id");
        if (p_ptr == NULL) {
            rebx_set_param_int(sim->extras,  (struct rebx_node**) &(p->ap), "fc_id", id_count);
            id_count++;
        }
        else{
            continue;
        }
    }
    rebx_set_param_int(sim->extras, &collision_resolve->ap, "fc_id_max", id_count);
}

static void output_collision_to_file(const char* filename, double t, enum COLLISION_TYPE collision_type, int new_id, int parent1_id, int parent2_id, double new_mass, double parent1_initial_mass, double parent2_initial_mass, double new_radius, double parent1_initial_radius, double parent2_initial_radius, double v_impact, double theta_impact){
    FILE* of = fopen(filename, "a");
    fprintf(of, "%e,", t);
    fprintf(of, "%d,", collision_type);
    fprintf(of, "%d,", new_id);
    fprintf(of, "%d,", parent1_id);
    fprintf(of, "%d,", parent2_id);
    fprintf(of, "%e,", new_mass);
    fprintf(of, "%e,", parent1_initial_mass);
    fprintf(of, "%e,", parent2_initial_mass);
    fprintf(of, "%e,", new_radius);
    fprintf(of, "%e,", parent1_initial_radius);
    fprintf(of, "%e,", parent2_initial_radius);
    fprintf(of, "%e,", v_impact);
    fprintf(of, "%e", theta_impact);
    fprintf(of, "\n");
    fclose(of);
}

// Function to merge two particles
static enum REB_COLLISION_RESOLVE_OUTCOME merge(struct reb_simulation* const sim, struct rebx_collision_resolve* const collision_resolve, struct reb_collision c, double v_impact, double theta_impact, enum COLLISION_TYPE collision_type){
    struct reb_particle* pi = &(sim->particles[c.p1]); // First object in collision
    struct reb_particle* pj = &(sim->particles[c.p2]); // Second object in collison

    double parent_1_initial_mass = pi->m;
    double parent_2_initial_mass = pj->m;
    double parent_1_initial_radius = pi->r;
    double parent_2_initial_radius = pj->r;

    double invmass = 1.0/(pi->m + pj->m);

    // Merge by conserving mass, volume and momentum
    pi->vx = (pi->vx*pi->m + pj->vx*pj->m)*invmass;
    pi->vy = (pi->vy*pi->m + pj->vy*pj->m)*invmass;
    pi->vz = (pi->vz*pi->m + pj->vz*pj->m)*invmass;
    pi->x  = (pi->x*pi->m + pj->x*pj->m)*invmass;
    pi->y  = (pi->y*pi->m + pj->y*pj->m)*invmass;
    pi->z  = (pi->z*pi->m + pj->z*pj->m)*invmass;
    double new_mass = pi->m + pj->m; // For printing
    pi->m  = new_mass;
    double new_radius = cbrt(pi->r*pi->r*pi->r + pj->r*pj->r*pj->r);
    pi->r = new_radius;
    pi->last_collision = sim->t;

    const char* particle_list_file = rebx_get_param(sim->extras, collision_resolve->ap, "fc_particle_list_file");
    if (particle_list_file != NULL) { // REBX parameter set?
        set_fc_ids(sim, collision_resolve);
        int parent_1_id = *(int*) rebx_get_param(sim->extras, pi->ap, "fc_id");
        int parent_2_id = *(int*) rebx_get_param(sim->extras, pj->ap, "fc_id");

        rebx_fragmenting_collisions_set_new_id(sim, collision_resolve, pi);
        int new_id = *(int*)rebx_get_param(sim->extras, pi->ap, "fc_id");
        output_collision_to_file(particle_list_file, sim->t, collision_type, new_id, parent_1_id, parent_2_id, new_mass, parent_1_initial_mass, parent_2_initial_mass, new_radius, parent_1_initial_radius, parent_2_initial_radius, v_impact, theta_impact); 
    }

    return REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P2; // Remove 2 particle from simulation
}

// Function to make fragments
static enum REB_COLLISION_RESOLVE_OUTCOME make_fragments(struct reb_simulation* const sim, struct rebx_collision_resolve* const collision_resolve,struct reb_collision c, double lr_mass, double slr_mass, double v_impact, double theta_impact, enum COLLISION_TYPE collision_type){
    // Get minimum fragment mass value
    // This is defined by the user in their setup
    double min_frag_mass;
    const double* min_frag_mass_ptr = rebx_get_param(sim->extras, collision_resolve->ap, "fc_min_frag_mass");
    if (min_frag_mass_ptr != NULL) {
        min_frag_mass = *min_frag_mass_ptr;
        // If it's valid, check if the value is not 0
        if (*min_frag_mass_ptr <= 0.0) {
            reb_simulation_error(sim, "Minimum fragment mass invalid (<= 0).\n");
            return REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
        }
    }else{
        reb_simulation_error(sim, "User needs to specify minimum fragment mass `fc_min_frag_mass`.\n");
        return REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
    } 
    double separation_distance_scale = 4; // Default value
    const double* separation_distance_scale_ptr = rebx_get_param(sim->extras, collision_resolve->ap, "fc_separation_distance_scale");
    if (separation_distance_scale_ptr != NULL) {
        separation_distance_scale = *separation_distance_scale_ptr; 
    } 

    struct reb_particle* pi = &(sim->particles[c.p1]); // First object in collision
    struct reb_particle* pj = &(sim->particles[c.p2]); // Second object in collison

    // Object with the higher mass will be the target, and object with lower mass will be the projectile
    struct reb_particle* target;     
    struct reb_particle* projectile; 

    enum REB_COLLISION_RESOLVE_OUTCOME outcome = REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
    if (pi->m >= pj->m){
        target = pi;    
        projectile = pj;
        outcome = REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P2;
    }else{
        target = pj;   
        projectile = pi;
        outcome = REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P1; 
    }

    struct reb_particle com = reb_particle_com_of_pair(*target, *projectile); // Center of mass (COM) of target and projectile

    double target_initial_mass = target->m; // Will use later for printing
    double projectile_initial_mass = projectile->m; // Will use later for printing
    double initial_mass = target->m + projectile->m; // initial mass of two colliders

    double target_initial_radius = target->r;
    double projectile_initial_radius = projectile->r;
    double r_tot = target->r + projectile->r; // Sum of radii or two colliders
    double remaining_mass = initial_mass - lr_mass - slr_mass; // Remaning mass, will turn into fragments
    double rho = target->m/(4./3*M_PI*pow(target ->r, 3)); // Target's density

    // slr_mass is the mass of the second largest remnant (refer to documentation for more info)
    // If slr_mass is non-zero, then we have a big fragment with mass slr_mass,
    // And a few small fragments.
    double n_big_frag = 0;
    if(slr_mass > 0){
        n_big_frag = 1;
    }

    /*
     * COMPUTING MASS OF FRAGMENTS
     */
    // We draw fragment masses from a power law based on LS2012
    double max_frag_mass = 0.5 * lr_mass;
    if(slr_mass > 0){
        max_frag_mass = 0.5 * slr_mass;
    }
    double powerlaw_slope = 3; // Arbitrary, from LS2012 table 1
    double m_frags_array[10000] = {0.0};
    double sum_m_frags = 0;
    int index = 0;
    int len_m_frags_array = 0;
    double ratio = 0;

    // Draw until we reach above the remaining mass. Then discard the last one and distribute its mass among others.
    while(sum_m_frags < remaining_mass && index < 10000){
        m_frags_array[index] = reb_random_powerlaw(sim, min_frag_mass, max_frag_mass, powerlaw_slope);
        sum_m_frags += m_frags_array[index];
        index += 1;
    }
    if(index >= 10000){
        reb_simulation_error(sim, "Number of fragments produced is above permitted value. Increase minimum fragment mass.\n");
        return REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
    }else if (index == 1){
        m_frags_array[0] = remaining_mass;
    }else if(index > 1){
        // Discard last fragment, and distribute the remaining mass between other fragments
        len_m_frags_array = index - 2; // These will be indexes of the frags we want to keep
        sum_m_frags = sum_m_frags - m_frags_array[len_m_frags_array + 1]; // Discard last fragment
        m_frags_array[len_m_frags_array + 1] = 0; 
        ratio = remaining_mass/sum_m_frags;
        // Distribute the mass between others
        for(int i=0; i<=len_m_frags_array; i++){
            m_frags_array[i] *= ratio;
        }
    }


    // n_frag is total number of fragments
    double n_frag = (len_m_frags_array + 1) + n_big_frag;

    /*
     * UPDATE TARGET TO BE THE LARGEST REMNANT
     */
    // We replace target with the largest remnant, and assign it the position and velocity of COM
    target -> last_collision = sim->t; // Update time of last collision
    target -> m = lr_mass; // Update target mass with lr_mass
    double lr_radius = get_radii(lr_mass, rho);
    target -> r = lr_radius; // Update target radius, keeping density constant
                             // Update target position with COM
    target->x = com.x; 
    target->y = com.y;
    target->z = com.z;
    // Update target velocity with COM velocity
    target->vx = com.vx;
    target->vy = com.vy;
    target->vz = com.vz;
    // Magnitude of lr_mass velocity, later to be used in computing fragment velocities
    double v_lr = sqrt((com.vx * com.vx) + (com.vy * com.vy) + (com.vz * com.vz));

    // Save parents IDs, to be printed later
    set_fc_ids(sim, collision_resolve);
    int parent_t_id = *(int*) rebx_get_param(sim->extras, target->ap, "fc_id");
    int parent_p_id = *(int*) rebx_get_param(sim->extras, projectile->ap, "fc_id");
    // Save new ID for lr
    rebx_fragmenting_collisions_set_new_id(sim, collision_resolve, target);
    int new_id = *(int*)rebx_get_param(sim->extras, target->ap, "fc_id");
    const char* particle_list_file = rebx_get_param(sim->extras, collision_resolve->ap, "fc_particle_list_file");
    if (particle_list_file != NULL) { // REBX parameter set?
        output_collision_to_file(particle_list_file, sim->t, collision_type, new_id, parent_t_id, parent_p_id, lr_mass, target_initial_mass, projectile_initial_mass, lr_radius, target_initial_radius, projectile_initial_radius, v_impact, theta_impact); 
    }

    //Define mxsum variable to keep track of center of mass (mass times position)
    struct reb_vec3d mxsum = {.x = 0, .y = 0, .z = 0};
    mxsum.x += target->m * target->x;
    mxsum.y += target->m * target->y;
    mxsum.z += target->m * target->z;

    //Define mvsum variable to keep track of momentum (mass times velocity)
    struct reb_vec3d mvsum = {.x = 0, .y = 0, .z = 0};
    mvsum.x += target->m * target->vx;
    mvsum.y += target->m * target->vy;
    mvsum.z += target->m * target->vz;

    /**
     * LOCATING FRAGMENTS 
     * Now we need to position fragments. Following Chambers (2013),
     * we first find a "collision plane" which is the plane crossing the two vectors of relative velocity
     * and relative position of target and projectile. 
     * Then, we draw a circle around the center of mass of target and projectile in this plane.
     * We position fragments in this circle, with equal angular seperation.
     * We assign them velocities relative to their masses, based on equipartition of kintetic energy. 
     * The velocity vector directions are derived following radii of the fragment circle, getting away from the COM.
     * For more information, refer to Chambers (2013) and Childs and Steffen (2022).
     */

    //Relative velocity between target and projectile in x, y, z
    struct reb_vec3d dv = {.x = target->vx - projectile->vx, .y = target->vy - projectile->vy, .z = target->vz - projectile->vz};

    //Relative location of target and projectile in x, y, z
    struct reb_vec3d dr = {.x = target->x - projectile->x, .y = target->y - projectile->y, .z = target->z - projectile->z};

    //Vector in the direction of relative velocity between target and projectile (unit vector of relative v)
    struct reb_vec3d unit_dv = reb_vec3d_normalize(dv);

    // normal_coll_plane : Vector normal to the collision plane
    //(cross product of rel. velocity and position between target and projectile)
    struct reb_vec3d normal_coll_plane = reb_vec3d_cross(dv, dr);

    //Turn the normal_coll_plane vector into a unit vector
    double normal_coll_plane_mag = sqrt(reb_vec3d_length_squared(normal_coll_plane));

    // Handle the case where the normal vector is zero (relative velocity and position are collinear)
    if (normal_coll_plane_mag < 1e-15) { 
        // Choose an arbitrary vector perpendicular to unit_dv
        struct reb_vec3d rand_vec_1 = {.x = 1.0, .y = 0.0, .z = 0.0};
        struct reb_vec3d rand_vec_2 = {.x = 0.0, .y = 1.0, .z = 0.0};
        struct reb_vec3d arb_vec = reb_vec3d_cross(unit_dv, rand_vec_1);

        if (reb_vec3d_length_squared(arb_vec) < 1e-15) { 
            arb_vec = reb_vec3d_cross(unit_dv, rand_vec_2);
        }
        normal_coll_plane = arb_vec;
        normal_coll_plane = reb_vec3d_normalize(normal_coll_plane);
    }
    // Normalize normal_coll_plane (now guaranteed non-zero)
    normal_coll_plane = reb_vec3d_normalize(normal_coll_plane);

    // normal_to_vrel : Vector normal to the relative velocity of target and projectile, in the collision plane
    // (cross product of normal_coll_plane and relative velocity)
    struct reb_vec3d normal_to_vrel = reb_vec3d_cross(dv, normal_coll_plane);
    normal_to_vrel = reb_vec3d_normalize(normal_to_vrel);

    // Separation distance is the distance between largest remnant and each fragment
    double separation_distance = separation_distance_scale * r_tot;

    // Seperation angle between fragments
    double theta_sep = (2.*M_PI)/n_frag;

    /*
     * ADD FRAGMENTS TO THE SIMULATION
     */
    // fragments are placed in the collision plane, in a circle with radius of separation distance.
    // Relative velocity unit vector and the vector orthogonal to that (normal_to_vrel) are used as
    // the reference frame to place fragments.
    // Fragments are placed with equal angular distances of each other (theta_sep).  

    // Add big fragment, if exists
    if (n_big_frag == 1){
        struct reb_particle big_frag = {0};
        big_frag.m = slr_mass;

        big_frag.x = com.x + separation_distance * unit_dv.x;
        big_frag.y = com.y + separation_distance * unit_dv.y;
        big_frag.z = com.z + separation_distance * unit_dv.z;

        double slr_v = v_lr * sqrt(lr_mass/slr_mass);

        big_frag.vx = com.vx + slr_v * unit_dv.x;
        big_frag.vy = com.vy + slr_v * unit_dv.y;
        big_frag.vz = com.vz + slr_v * unit_dv.z;

        double slr_radius = get_radii(slr_mass, rho);
        big_frag.r = slr_radius;

        // Record collision
        big_frag.last_collision = sim->t;

        // Add to mxsum vector to keep track of COM
        mxsum.x +=big_frag.m * big_frag.x;
        mxsum.y += big_frag.m * big_frag.y;    
        mxsum.z += big_frag.m * big_frag.z;

        // Add to mvsum vector to keep track of momentum
        mvsum.x += big_frag.m * big_frag.vx;
        mvsum.y += big_frag.m * big_frag.vy;    
        mvsum.z += big_frag.m * big_frag.vz;

        // Add particle to simulation
        reb_simulation_add(sim, big_frag);

        // Save new ID with parents to particle ID list
        rebx_fragmenting_collisions_set_new_id(sim, collision_resolve, &sim->particles[sim->N - 1]);
        const char* particle_list_file = rebx_get_param(sim->extras, collision_resolve->ap, "fc_particle_list_file");
        if (particle_list_file != NULL) { // REBX parameter set?
            struct reb_particle* newly_added_particle = &(sim->particles[sim->N - 1]); 
            int new_id = *(int*)rebx_get_param(sim->extras, newly_added_particle->ap, "fc_id");
            output_collision_to_file(particle_list_file, sim->t, collision_type, new_id, parent_t_id, parent_p_id, slr_mass, target_initial_mass, projectile_initial_mass, slr_radius, target_initial_radius, projectile_initial_radius, v_impact, theta_impact); 
        }
    }

    // Add small fragments
    // j = 0 is reserved for big_frag, if exists
    for (int j=1; j <= n_frag - n_big_frag; j++){          
        struct reb_particle fragment = {0};
        fragment.m = m_frags_array[j-1]; 

        fragment.x = com.x + separation_distance*(cos(theta_sep*j)*unit_dv.x + sin(theta_sep*j)*normal_to_vrel.x);
        fragment.y = com.y + separation_distance*(cos(theta_sep*j)*unit_dv.y + sin(theta_sep*j)*normal_to_vrel.y);
        fragment.z = com.z + separation_distance*(cos(theta_sep*j)*unit_dv.z + sin(theta_sep*j)*normal_to_vrel.z);

        double frag_velocity = v_lr * sqrt(lr_mass/fragment.m);

        fragment.vx = com.vx + frag_velocity*(cos(theta_sep*j)*unit_dv.x + sin(theta_sep*j)*normal_to_vrel.x);
        fragment.vy = com.vy + frag_velocity*(cos(theta_sep*j)*unit_dv.y + sin(theta_sep*j)*normal_to_vrel.y);
        fragment.vz = com.vz + frag_velocity*(cos(theta_sep*j)*unit_dv.z + sin(theta_sep*j)*normal_to_vrel.z);

        // Fragment radius is derived based on target's density
        double targ_rho = target->m/(4./3*M_PI*pow(target->r,3));
        fragment.r = get_radii(m_frags_array[j-1], targ_rho);

        // Record collision
        fragment.last_collision = sim->t;


        // Add to mxsum vector to keep track of COM
        mxsum.x +=fragment.m*fragment.x;
        mxsum.y += fragment.m*fragment.y;    
        mxsum.z += fragment.m*fragment.z;

        // Add to mvsum vector to keep track of momentum
        mvsum.x += fragment.m*fragment.vx;
        mvsum.y += fragment.m*fragment.vy;    
        mvsum.z += fragment.m*fragment.vz;

        // Finally add fragment to simulation.
        reb_simulation_add(sim, fragment); 

        // Set the fragment ID
        rebx_fragmenting_collisions_set_new_id(sim, collision_resolve, &sim->particles[sim->N - 1]);
        struct reb_particle* newly_added_particle = &(sim->particles[sim->N - 1]); //First object in collision
        int new_id = *(int*)rebx_get_param(sim->extras, newly_added_particle->ap, "fc_id");

        // Save fragment ID into particle ID list
        const char* particle_list_file = rebx_get_param(sim->extras, collision_resolve->ap, "fc_particle_list_file");
        if (particle_list_file != NULL) { // REBX parameter set?
            output_collision_to_file(particle_list_file, sim->t, collision_type, new_id, parent_t_id, parent_p_id, fragment.m, target_initial_mass, projectile_initial_mass, fragment.r, target_initial_radius, projectile_initial_radius, v_impact, theta_impact); 
        }

    }

    // Now we correct for the COM and momentum offsets.
    // First, we need to see how much we are off from the COM in the begining. 
    struct reb_vec3d xoff = {.x = com.x - mxsum.x/initial_mass, .y = com.y - mxsum.y/initial_mass, .z = com.z - mxsum.z/initial_mass};
    struct reb_vec3d voff = {.x = com.vx - mvsum.z/initial_mass, .y = com.vy - mvsum.y/initial_mass, .z = com.vz - mvsum.z/initial_mass};

    // Reassign new position and velocity to target (who is replaced by the largest remnant) to account for offsets:
    target -> x +=  xoff.x*target->m/initial_mass;
    target -> y += xoff.y*target->m/initial_mass; 
    target -> z += xoff.z*target->m/initial_mass; 
    target -> vx += voff.x*target->m/initial_mass; 
    target -> vy += voff.y*target->m/initial_mass; 
    target -> vz += voff.z*target->m/initial_mass; 

    // Reassign position and velocity to fragments to correct for offsets.
    for (int i = sim->N - n_frag; i < sim->N; i++){ 
        // mass fraction of fragment versus total initial mass
        double mass_fraction = sim->particles[i].m/initial_mass;

        sim->particles[i].x += xoff.x*mass_fraction;
        sim->particles[i].y += xoff.y*mass_fraction;
        sim->particles[i].z += xoff.z*mass_fraction;

        sim->particles[i].vx += voff.x*mass_fraction;
        sim->particles[i].vy += voff.y*mass_fraction;
        sim->particles[i].vz += voff.z*mass_fraction;
    }

    return outcome;
}

/*
 * Main function to decide the collision outcome, derive new masses, positions and velocities.
 * Equations are derived from LS2012 and Chambers (2013).
 */
enum REB_COLLISION_RESOLVE_OUTCOME rebx_fragmenting_collisions(struct reb_simulation* const sim, struct rebx_collision_resolve* const collision_resolve, struct reb_collision c){
    // Setting minimum fragment mass
    const double* min_frag_mass_ptr = rebx_get_param(sim->extras, collision_resolve->ap, "fc_min_frag_mass");
    double min_frag_mass;
    if (min_frag_mass_ptr != NULL) {
        // If it's valid, check if the value is not 0
        min_frag_mass = *min_frag_mass_ptr;
        if (*min_frag_mass_ptr <= 0.0) {
            reb_simulation_error(sim, "Minimum fragment mass invalid (<= 0).\n");
            return REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
        }
    }else{
        reb_simulation_error(sim, "User needs to specify minimum fragment mass `fc_min_frag_mass`.\n");
        return REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
    }   
    double rho1 = 1.684e6; // Default value. Units of Msun/AU^3 
    const double* rho1_ptr = rebx_get_param(sim->extras, collision_resolve->ap, "fc_rho1");
    if (rho1_ptr != NULL) {
        rho1 = *rho1_ptr; 
    } 
    double cstar = 1.8; // Default value 
    const double* cstar_ptr = rebx_get_param(sim->extras, collision_resolve->ap, "fc_cstar");
    if (cstar_ptr != NULL) {
        cstar = *cstar_ptr; 
    } 

    struct reb_particle* pi = &(sim->particles[c.p1]); // First object in collision
    struct reb_particle* pj = &(sim->particles[c.p2]); // Second object in collison

    if (pi->last_collision==sim->t || pj->last_collision==sim->t){
        return 0;
    }

    // Object with the higher mass will be the target, and object with lower mass will be the projectile
    struct reb_particle* target;     
    struct reb_particle* projectile; 

    // Object with the higher mass will be the target, and object with lower mass will be the projectile
    if (pi->m >= pj->m){
        target = pi;    
        projectile = pj;
    }else{
        target = pj;   
        projectile = pi; 
    }

    if (target->m == 0.0){
        reb_simulation_error(sim, "Target mass is zero.\n");
        return REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
    }
    if (projectile->m == 0.0){
        reb_simulation_error(sim, "Projectile mass is zero.\n");
        return REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
    }

    // Some useful parameters
    const double target_initial_mass = target->m; // To be printed 
    const double projectile_initial_mass = projectile->m; // To be printed 
    const double initial_mass = target_initial_mass + projectile_initial_mass; // Initial mass of two colliders
    const double target_initial_radius = target->r;
    const double projectile_initial_radius = projectile->r;
    const double r_tot = target_initial_radius + projectile_initial_radius; // Sum of radii of colliders

    // Relative positions
    struct reb_vec3d r = {.x = target->x - projectile->x, .y = target->y - projectile->y, .z = target->z - projectile->z};
    double distance_mag = sqrt(reb_vec3d_length_squared(r));

    // Relative velocities
    struct reb_vec3d dv = {.x = target->vx - projectile->vx, .y = target->vy - projectile->vy, .z = target->vz - projectile->vz};
    double dv_mag = sqrt(reb_vec3d_length_squared(dv));

    // Angular momentum vector, cross product of relative velocity and position
    struct reb_vec3d h = reb_vec3d_cross(dv, r);
    double h_mag = sqrt(reb_vec3d_length_squared(h));

    // Impact velocity (refer to eq. 1 in Childs and Steffen (2022))
    double v_imp = sqrt(dv_mag*dv_mag + 2 * (sim->G) * initial_mass * (1./r_tot - 1./distance_mag));

    // If collision is detected after physical contact,
    // then the distance between two objects is less than sum of their radii. 
    // In this case impact velocity is just relative velocity.
    if (1./r_tot < 1./distance_mag){
        v_imp = dv_mag;
    }

    // Impact parameter, defined as b = (R_t + R_p)sin(theta). Refer to Figure 2. in LS2012.
    double b = h_mag/v_imp;
    double theta_i = asin(b/distance_mag); // Impact angle (radians), to be printed later

    if (isnan(b)){
        reb_simulation_error(sim, "b is not a number.");
        return REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
    }

    // The following are steps to find collision energy, and derive largest remnant mass accordingly
    // Refer to LS2012 for the full description.
    // Refer to Chambers (2013) for a shortened description of equations. 
    // Chambers (2013) Eq. 2, reduced mass
    double mu = (target->m * projectile->m)/initial_mass;

    // LS2012 Eq. 7, the projected length of the projectile overlapping the target
    double l = r_tot-b;
    l = MIN(l, 2*projectile->r);

    // LS2012 Eq. 11, interacting mass fraction
    double alpha = (pow(l,2)*(3*projectile->r - l))/(4*pow(projectile->r, 3));
    alpha = MIN(1., alpha);

    // Specific energy per unit mass Q, Chambers (2013) eq. 1.
    double Q = 0.5 * pow(v_imp,2) * target->m * projectile->m / pow(initial_mass,2);

    // Mutual escape velocity of target and projectile
    double v_esc = pow(2.*(sim->G)*initial_mass/r_tot, 0.5);

    // LS2012 Eq. 12, reduced interacting mass for oblique impacts.
    double alphamu = (alpha * target->m * projectile->m)/(alpha * projectile->m + target->m);

    // Mass ratio of target and projectile, Chambers (2013) eq. 6
    double gamma = projectile->m/target->m;  

    // Chambers (2013) Eq. 4, combined radius of target and projectile with constant density
    double Rc1 = pow((initial_mass * 3)/(4. * M_PI * rho1), 1./3.);  

    // Chambers (2013) Eq. 3, critical value of impact energy for head-on collisions
    double Q0 = 0.8 * cstar * M_PI * rho1 * (sim->G) * pow(Rc1,2); 

    // Chambers (2013) Eq. 5, critical impact energy for oblique or different mass collisons.  
    double Q_star = pow(mu/alphamu, 1.5)*(pow(1+gamma, 2)/ (4*gamma))*Q0; 
    if (alpha == 0.0){
        reb_simulation_error(sim, "alpha (interacting mass fraction) = 0");
        return REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
    }

    // For equal mass and head-on collisions Q* = Q0.
    if (b == 0 && target->m == projectile->m){
        Q_star = Q0;
    }
    // Mass of largest remnant is derived based on Q and Q* ratio. (Chambers (2013) eq. 8)
    double lr_mass;
    double qratio = Q/Q_star;
    if (qratio < 1.8){
        lr_mass = initial_mass*(1.0-.5*qratio);
    }else{
        lr_mass = 0.1 * initial_mass * pow(qratio/1.8, -1.5);  
    }

    enum COLLISION_TYPE collision_type;
    enum REB_COLLISION_RESOLVE_OUTCOME outcome = REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;

    /*
     * DECIDE WHAT TO DO AFTER THE COLLISION
     */
    // Refer to documentation for the decision tree flowchart and description. 

    // If v_imp <= v_esc, merge.
    if (v_imp <= v_esc){
        collision_type = COLLISION_TYPE_MERGE;
        outcome = merge(sim, collision_resolve, c, v_imp, theta_i, collision_type);
        printf("Merging collision detected. (Case A)\n");
    }else{
        if(b >= target->r){ // Grazing regime
                            // Target's density
            double targ_rho = target->m/(4./3*M_PI*pow(target->r,3));

            // phi helps with finding part of the projectile that is NOT crossing the target
            double phi = 2*acos((l-projectile->r)/projectile->r);

            // LS2012 Eq. 46; cross section of projectile interacting with the target
            double A_interact = pow(projectile->r, 2)*((M_PI-(phi-sin(phi))/2.));  

            // LS2012 Eq. 47, interacting length
            double L_interact = 2.*pow(pow(target->r,2)-(pow(target->r-l/2.,2)), .5);

            // LS2012 Eq. 48, used in Chambers Eq. 11
            double beta = ((A_interact*L_interact) * targ_rho)/target->m;

            // Based on Chambers Eq. 11, subscript g refers to "grazing"
            double Rc1_g = pow(3./(4.*M_PI*rho1)*(beta * target->m + projectile->m), 1./3.);

            // Chambers Eq. 11
            double Q0_g = .8*cstar*M_PI*rho1*sim->G*pow(Rc1_g, 2); 

            // LS2012 Eq. 46-59
            double gamma_g = beta * target->m/projectile->m;

            // Chambers Eq. 10
            double Q_star_g = (pow(1+gamma_g, 2)/4*gamma_g)* Q0_g; 

            // Chambers Eq. 13
            double mu_g = (beta * target->m*projectile->m)/(beta * target->m+projectile->m);  

            // Chambers Eq. 12
            double Q_g = .5*(mu_g*pow(v_imp,2))/(beta*target->m+projectile->m); 

            /* If  velocity in the hit-and-run regime is very low, the collision
             * might eventually lead to a merger. Here, we compute the threshhold velocity for this event,
             * called critical velocity. If v < v_crit, then we have a "graze and merge" event.
             */

            // c1 to c4 are constants used in Chambers Eq. 17
            double c1 = 2.43; 
            double c2 = -0.0408;
            double c3 = 1.86;
            double c4 = 1.08;

            // Chambers eq. 16
            double zeta = pow((1 - gamma)/(1 + gamma),2);

            // This helps with writing Chambers eq. 15
            double fac = pow(1-b/(target->r + projectile->r),2.5);

            // Velocity threshhold between graze-and-merge and hit-and-run, Chambers Eq. 15
            double v_crit = v_esc*(c1*zeta*fac + c2*zeta +c3*fac + c4);

            // If impact velocity is less than v_crit, we have graze-and-merge
            if (v_imp <= v_crit){       
                collision_type = COLLISION_TYPE_MERGE_B;
                printf("Merging collision detected. (Case B)\n");
                outcome = merge(sim, collision_resolve, c, v_imp, theta_i, collision_type);
            }else{
                if(lr_mass < target->m){
                    if((target->m + projectile->m - lr_mass) < min_frag_mass){
                        collision_type = COLLISION_TYPE_BOUNCE_F;
                        printf("Elastic bounce, (Case F).\n");
                        outcome = REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
                        reb_collision_resolve_hardsphere(sim,c);
                    }else{
                        collision_type = COLLISION_TYPE_GRAZING_G;
                        printf("Grazing erosion, (Case G).\n");
                        outcome = make_fragments(sim, collision_resolve, c, lr_mass, 0, v_imp, theta_i, collision_type);
                    }
                }else{
                    // lr_dag_mass : Second largest remnant mass, Chambers Eq. 14
                    double lr_dag_mass;
                    if (Q_g < 1.8*Q_star_g){
                        lr_dag_mass = (beta*target->m + projectile->m)*(1 - Q_g/ (2*Q_star_g));
                    }else{
                        lr_dag_mass = (beta*target->m + projectile->m)/10 * pow(Q_g/(1.8*Q_star_g), -1.5);
                    }
                    if(lr_dag_mass < min_frag_mass){
                        collision_type = COLLISION_TYPE_BOUNCE_H;
                        printf("Elastic bounce, (Case H).\n");
                        outcome = REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
                        reb_collision_resolve_hardsphere(sim,c);
                    }else{
                        if((target->m + projectile->m - lr_mass - lr_dag_mass) < min_frag_mass){
                            collision_type = COLLISION_TYPE_BOUNCE_I;
                            printf("Elastic bounce, (Case I).\n");
                            outcome = REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE;
                            reb_collision_resolve_hardsphere(sim,c);
                        }else{
                            collision_type = COLLISION_TYPE_HITANDRUN;
                            printf("Hit-and-run, (Case J).\n");
                            outcome = make_fragments(sim, collision_resolve, c, lr_mass, lr_dag_mass, v_imp, theta_i, collision_type);
                        }
                    }
                }
            }
        } 
        else{ //non-grazing regime
            if (initial_mass - lr_mass < min_frag_mass){ //Not meeting minimum fragment mass threshold
                collision_type = COLLISION_TYPE_MERGE_C;
                printf("Non grazing, M_rem to small. Merging collision detected. (Case C)\n");
                outcome = merge(sim, collision_resolve, c, v_imp, theta_i, collision_type);
            }else{ //Can make fragments. lr_mass can be larger or smaller than the target
                if(lr_mass > target->m){
                    collision_type = COLLISION_TYPE_ACCRETION;
                    printf("Non grazing, lr_mass > M_t. Accretion. (Case 4, D)\n");
                    outcome = make_fragments(sim, collision_resolve, c, lr_mass, 0, v_imp, theta_i, collision_type);
                }else{
                    if(lr_mass < min_frag_mass){
                        collision_type = COLLISION_TYPE_SUPERCATASTROPHIC;
                        lr_mass = min_frag_mass;
                        printf("Super Catastrophic, lr_mass < min_frag_mass (Case E)\n");
                        outcome = make_fragments(sim, collision_resolve, c, lr_mass, 0, v_imp, theta_i, collision_type);
                    }else{
                        collision_type = COLLISION_TYPE_ACCRETION;
                        printf("Non grazing, lr_mass < M_t. Erosion. (Case 4, D)\n");
                        outcome = make_fragments(sim, collision_resolve, c, lr_mass, 0, v_imp, theta_i, collision_type);
                    }
                }
            }
        }
    }
    // Print collision data for elastic bounces
    if (collision_type == COLLISION_TYPE_BOUNCE_F || collision_type == COLLISION_TYPE_BOUNCE_H || collision_type == COLLISION_TYPE_BOUNCE_I){
        const char* particle_list_file = rebx_get_param(sim->extras, collision_resolve->ap, "fc_particle_list_file");
        if (particle_list_file != NULL) { // REBX parameter set?
            set_fc_ids(sim, collision_resolve);

            int parent_t_id = *(int*)rebx_get_param(sim->extras, target->ap, "fc_id");
            int parent_p_id = *(int*)rebx_get_param(sim->extras, projectile->ap, "fc_id");
            // Target and projectile both get new IDs
            rebx_fragmenting_collisions_set_new_id(sim, collision_resolve, target);
            rebx_fragmenting_collisions_set_new_id(sim, collision_resolve, projectile);
            int new_id_t = *(int*)rebx_get_param(sim->extras, target->ap, "fc_id");
            int new_id_p = *(int*)rebx_get_param(sim->extras, projectile->ap, "fc_id");

            // Print data for object 1 (target)
            output_collision_to_file(particle_list_file, sim->t, collision_type, new_id_t, parent_t_id, parent_p_id, target_initial_mass, target_initial_mass, projectile_initial_mass, target_initial_radius, target_initial_radius, projectile_initial_radius, v_imp, theta_i); 

            // Print data for object 2 (projectile)
            output_collision_to_file(particle_list_file, sim->t, collision_type, new_id_p, parent_t_id, parent_p_id, projectile_initial_mass, target_initial_mass, projectile_initial_mass, projectile_initial_radius, target_initial_radius, projectile_initial_radius, v_imp, theta_i); 
        }
    }

    return outcome;
} 
