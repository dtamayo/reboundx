/**
 * @file    merging_collisions.c
 * @brief   Adds a REBOUNDx collision module which always merges particles
 * @author  Hanno Rein <hanno.rein@utoronto.ca>
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
 * $Fragmenting collisions$     // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script). 
 * 
 * ======================= ===============================================
 * Authors                 H. Rein 
 * Based on                None
 * C Example               :ref:`c_fragmenting_collisions`
 * ======================= ===============================================
 * 
 * This is a simple example implementation of a REBOUNDx collision module.
 * The outcome is a largest remnant and multiple fragments.
 * 
 * **Effect Parameters**
 * 
 * *None*
 * 
 * **Particle Parameters**
 * 
 * *None*
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

/** 
* Function to get cross product of two vectors
* First vector is (Ax, Ay, Ax) and second one is (Bx, By, Bz)
* Results will be saved in resultX, resultY and resultZ 
*/
void get_cross_product(double Ax, double Ay, double Az,
                             double Bx, double By, double Bz,
                             double* resultX, double* resultY, double* resultZ) {
    *resultX = Ay * Bz - Az * By;
    *resultY = Az * Bx - Ax * Bz;
    *resultZ = Ax * By - Ay * Bx;
}

/**
 * Function to get magnitude of a vector.
 */
double get_mag(double x, double y, double z){
    return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
}  


int rebx_fragmenting_collisions(struct reb_simulation* const sim, struct rebx_collision_resolve* const collision_resolve, struct reb_collision c){
    struct reb_particle* pi = &(sim->particles[c.p1]); //First object in collision
    struct reb_particle* pj = &(sim->particles[c.p2]); //Second object in collison

    //Object with the higher mass will be the target, and object with lower mass will be the projectile
    struct reb_particle* target;     
    struct reb_particle* projectile; 

    //Object with the higher mass will be the target, and object with lower mass will be the projectile
    if (pi->m >= pj->m){
        target = pi;    
        projectile = pj;
    }
    else{
        target = pj;   
        projectile = pi; 
    }

    struct reb_particle com = reb_particle_com_of_pair(*target, *projectile); //Center of mass (COM) of target and projectile
    double initial_mass = target -> m + projectile -> m; //Sum of mass of two objects
    double Mlr; //Mass of the largest remnant
    double remaining_mass = initial_mass - Mlr; //Remaning mass, will turn into fragments
    double rho = target->m/(4./3*M_PI*pow(target ->r, 3)); //Target's density
    double min_frag_mass;

    //NOTE: There is a note about hit-and-run collisions and tracking the "big fragment" but we ignore this for now
    
    //We divide the remaning mass into equal mass fragments
    int n_frag = remaining_mass/min_frag_mass; //number of fragments
    double m_frag = remaining_mass/n_frag; //mass of each fragment

    //Define mxsum variable to keep track of center of mass (mass times position)
    double mxsum[3] = {0, 0, 0}; // For x, y, z

    //Define mvsum variable to keep track of momentum (mass times velocity)
    double mvsum[3] = {0, 0, 0}; // For x, y, z

    //We replace target with the largest remnant, and assign it the position and velocity of COM
    target -> last_collision = sim->t; //Update time of last collision
    target -> m = Mlr; //Update target mass with Mlr
    target -> r = get_radii(Mlr, rho); //Update target radius, keeping density
    //Update target position with COM
    target->x = com.x; 
    target->y = com.y;
    target->z = com.z;
    //Update target velocity with COM velocity
    target->vx = com.vx;
    target->vy = com.vy;
    target->vz = com.vz;

    //NOTE: Childs code swaps Mlr and fragment, if Mlr falls below min_frag_mass. Need to address this later.

    //Track center of mass    
    mxsum[0] += target->m * target->x;
    mxsum[1] += target->m * target->y;
    mxsum[2] += target->m * target->z;
    
    //Track momentum
    mvsum[0] += target->m * target->vx;
    mvsum[1] += target->m * target->vy;
    mvsum[2] += target->m * target->vz;

    /**
    * Now we need to position fragments. Following Chambers (2013),
    * we first find a "collision plane" which is the plane crossing the two vectors of relative velocity
    * and relative position of target and projectile. 
    * Then, we draw a circle around the center of mass of target and projectile in this plane.
    * We position fragments in this circle, with equal angular seperation.
    * We assign them velocities with magnitudes 5 percent higher than target and projectile's escape velocity. 
    * The velocity vectors are derived following radii of the fragment circle, getting away from the COM.
    * For more information, refer to Chambers (2013) and Childs and Steffen (2022).
    */

    //Seperation angle between fragments
    double theta_inc = (2.*M_PI)/n_frag;

    //Relative velocity between target and projectile in x, y, z
    double dvx = target->vx - projectile->vx;
    double dvy = target->vy - projectile->vy;
    double dvz = target->vz - projectile->vz;
    double dv_mag = get_mag(dvx, dvy, dvz);

    //Relative location of target and projectile in x, y, z
    double dx = target->x - projectile->x;
    double dy = target->y - projectile->y;
    double dz = target->z - projectile->z;
    double distance_mag = get_mag(dx, dy, dz);

    //Vector in the direction of relative velocity between target and projectile (unit vector of relative v)
    double unit_dvx = dvx/dv_mag;
    double unit_dvy = dvy/dv_mag;
    double unit_dvz = dvz/dv_mag;

    //Vector normal to the collision plane
    //(cross product of rel. velocity and position between target and projectile)
    double normal_coll_plane[3] = {0, 0, 0}; 
    get_cross_product(dvx, dvy, dvz, dx, dy, dz, &normal_coll_plane[0], &normal_coll_plane[1], &normal_coll_plane[2]);
    //Turn the normal_coll_plane vector into a unit vector
    double normal_coll_plane_mag = get_mag(normal_coll_plane[0], normal_coll_plane[1], normal_coll_plane[2]);
    normal_coll_plane[0] = normal_coll_plane[0]/normal_coll_plane_mag;
    normal_coll_plane[1] = normal_coll_plane[1]/normal_coll_plane_mag;
    normal_coll_plane[2] = normal_coll_plane[2]/normal_coll_plane_mag;


    //Vector normal to the relative velocity of target and projectile, in the collision plane
    //(cross product of normal_coll_plane and relative velocity)
    double normal_to_vrel[3] = {0, 0, 0};
    get_cross_product(dvx, dvy, dvz, normal_coll_plane[0], normal_coll_plane[1], normal_coll_plane[2],
                      &normal_to_vrel[0], &normal_to_vrel[1], &normal_to_vrel[2]);
    //Turn it into a unit vector
    double normal_to_vrel_mag = get_mag(normal_to_vrel[0], normal_to_vrel[1], normal_to_vrel[2]);
    normal_to_vrel[0] = normal_to_vrel[0]/normal_to_vrel_mag;
    normal_to_vrel[1] = normal_to_vrel[1]/normal_to_vrel_mag;
    normal_to_vrel[2] = normal_to_vrel[2]/normal_to_vrel_mag;

    //NOTE: continue with fragment velocity, then adding fragments


    return 2; // Remove 2 particle from simulation
}

