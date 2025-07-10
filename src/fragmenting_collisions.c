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
* Function to get cross product of two vectors.
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

/*
Function to get dot product of two vectors.
*/
double get_dot(double x1, double y1, double z1, double x2, double y2, double z2){ 
    return (x1*x2)+(y1*y2)+(z1*z2);
}

/**
 * Function to get magnitude of a vector.
 */
double get_mag(double x, double y, double z){
    return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
}  

/*
* Function to get radius of an object, given mass and density
*/
double get_radii(double m, double rho){
    return pow((3.*m)/(4.*M_PI*rho),1./3.);
}   

//Global parameters, need to be used defined
double separation_distance_scale = 4;
double min_frag_mass = 1.4e-8;

/*
* Function to get mass of the largest remnant, based on the energy involved in the collision.
* Equations are derived from Leinhardt and Stewart (2011) and Chambers (2013).
*/
double get_mass_of_largest_remnant(struct reb_simulation* const sim, struct reb_collision c){
    double Mlr;
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

    //Some useful parameters
    double initial_mass = target->m + projectile->m; //initial mass of two colliders
    double r_tot = target->r + projectile->r; //sum of radii of colliders
    double G = sim->G; //gravitational constant

    //Relative positions
    double dx = target->x - projectile->x;
    double dy = target->y - projectile->y;
    double dz = target->z - projectile->z;
    double distance_mag = get_mag(dx, dy, dz);

    //Relative velocities
    double dvx = target->vx - projectile->vx;
    double dvy = target->vy - projectile->vy;
    double dvz = target->vz - projectile->vz;
    double dv_mag = get_mag(dvx, dvy, dvz);

    //Angular momentum vector, cross product of relative velocity and position
    double hx, hy, hz;
    get_cross_product(dvx, dvy, dvz, dx, dy, dz, &hx, &hy, &hz);
    double h_mag = get_mag(hx, hy, hz);

    //Impact velocity (refer to eq. 1 in Childs and Steffen (2022))
    double v_imp = sqrt(dv_mag*dv_mag + 2 * G * initial_mass * (1./r_tot - 1./distance_mag));
    //If collision is detected after physical contact,
    //then the distance between two objects is less than sum of their radii. 
    //In this case impact velocity is just relative velocity.
    if (1./r_tot < 1./distance_mag){
        v_imp = dv_mag;
    }

    //Impact parameter. Refer to Figure 2. in Leinhardt and Stewart 2011.
    double b = h_mag/v_imp; 
    if (isnan(b)){
        reb_simulation_error(sim, "b is not a number.");
        return 0;
    }

    //The following are steps to find collision energy, and derive largest remnant mass accordingly
    //Refer to Leinhardt and Stewart (2011) for the full description.
    //Refer to Chambers (2013) for a shortened description of equations. 
    //Chambers (2013) Eq. 2, reduced mass
    double mu = (target->m * projectile->m)/initial_mass;

    //Leinhardt and Stewart (2011) Eq. 7, the projected length of the projectile overlapping the target
    double l = r_tot-b;
    l = MIN(l, 2*projectile->r);

    //Leinhardt and Stewart (2011) Eq. 11, interacting mass fraction
    double alpha = (pow(l,2)*(3*projectile->r - l))/(4*pow(projectile->r, 3));
    alpha = MIN(1., alpha);

    //Specific energy per unit mass Q, Chambers (2013) eq. 1.
    double Q = 0.5 * pow(v_imp,2) * target->m * projectile->m / pow(initial_mass,2);

    //Mutual escape velocity of target and projectcile
    double V_esc = pow(2.*G*initial_mass/r_tot, 0.5);

    //Leinhardt and Stewart (2011) Eq. 12, reduced interacting mass for oblique impacts.
    double alphamu = (alpha * target->m * projectile->m)/(alpha * projectile->m + target->m);

    //Mass ratio of target and projectile, Chambers (2013) eq. 6
    double gamma = projectile->m/target->m;  

    double Rc1 = pow((M_tot*3)/(4.*M_PI*rho1), 1./3.);  //Chambers Eq. 4, combined radius of target and projectile with constant density
    double Q0 = .8*cstar*M_PI*rho1*G*pow(Rc1,2);  //Chambers Eq. 3, critical value of impact energy for head-on collisions
    double Q_star = pow(mu/alphamu, 1.5)*(pow(1+gamma, 2)/ (4*gamma))*Q0;  //Chambers Eq. 5, critical value for oblique or different mass collisons.  
    if (alpha == 0.0){
        reb_simulation_error(r, "alpha = 0");
        return 0;
    }

    return Mlr;
}



int rebx_fragmenting_collisions(struct reb_simulation* const sim, struct rebx_collision_resolve* const collision_resolve, struct reb_collision c, double largest_remnant_mass){
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
    double initial_mass = target->m + projectile->m; //initial mass of two colliders
    double r_tot = target->r + projectile->r; //Sum of radii or two colliders
    double Mlr = largest_remnant_mass; //Mass of the largest remnant
    double remaining_mass = initial_mass - Mlr; //Remaning mass, will turn into fragments
    double rho = target->m/(4./3*M_PI*pow(target ->r, 3)); //Target's density

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

    // Handle the case where the normal vector is zero (relative velocity and position are collinear)
    if (normal_coll_plane_mag < 1e-15) { 
        // Choose an arbitrary vector perpendicular to unit_dv
        double arb_vec_x, arb_vec_y, arb_vec_z;
        get_cross_product(unit_dvx, unit_dvy, unit_dvz, 1.0, 0.0, 0.0, &arb_vec_x, &arb_vec_y, &arb_vec_z);
        if (get_mag(arb_vec_x, arb_vec_y, arb_vec_z) < 1e-15) { 
            get_cross_product(unit_dvx, unit_dvy, unit_dvz, 0.0, 1.0, 0.0, &arb_vec_x, &arb_vec_y, &arb_vec_z);
        }
        normal_coll_plane[0] = arb_vec_x;
        normal_coll_plane[1] = arb_vec_y;
        normal_coll_plane[2] = arb_vec_z;
        normal_coll_plane_mag = get_mag(normal_coll_plane[0], normal_coll_plane[1], normal_coll_plane[2]);
    }
    
    // Normalize normal_coll_plane (now guaranteed non-zero)
    normal_coll_plane[0] /= normal_coll_plane_mag;
    normal_coll_plane[1] /= normal_coll_plane_mag;
    normal_coll_plane[2] /= normal_coll_plane_mag;


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

    //Compute magnitude of fragment velocity. Here, we choose 5% more than the escape velocity.
    double G = sim->G;
    //escape velocity = (G.M_total/R_total)^(1/2)
    double V_esc = pow(2.*G*(initial_mass)/(r_tot), .5);
    //Separation distance is the distance between largest remnant and each fragment
    double separation_distance = separation_distance_scale * r_tot;
    //Fragment velocity, refer to Childs and Steffen (2022) eq. 1 for a similar computation of impact velocity. 
    //In summary, we need to subtract the potential energy, which is different at the moment of contact (where the 
    //distance between two particles is r_tot) and when the fragments are placed and leaving the largest remnant
    //(where the distance between lr and frag is = separation distance).
    double frag_velocity =sqrt(1.1*pow(V_esc,2) - 2 * G* initial_mass * (1./(r_tot) - 1./(separation_distance)));
    //Seperation angle between fragments
    double theta_sep = (2.*M_PI)/n_frag;

    //Add fragments to the simulation
    for (int j=1; j < n_frag + 1; j++){          
        struct reb_particle fragment = {0};

        //fragment mass computed above.
        fragment.m = m_frag; 
        //fragments are placed in the collision plane, in a circle with radius of separation distance.
        //Relative velocity unit vector and the vector orthogonal to that (normal_to_vrel) are used as
        //the reference frame to place fragments.
        //Fragments are placed with equal angular distances of each other (theta_sep).                
        fragment.x = com.x + separation_distance*(cos(theta_sep*j)*unit_dvx + sin(theta_sep*j)*normal_to_vrel[0]);
        fragment.y = com.y + separation_distance*(cos(theta_sep*j)*unit_dvy + sin(theta_sep*j)*normal_to_vrel[1]);
        fragment.z = com.z + separation_distance*(cos(theta_sep*j)*unit_dvz + sin(theta_sep*j)*normal_to_vrel[2]);
        fragment.vx = com.vx + frag_velocity*(cos(theta_sep*j)*unit_dvx + sin(theta_sep*j)*normal_to_vrel[0]);
        fragment.vy = com.vy + frag_velocity*(cos(theta_sep*j)*unit_dvy + sin(theta_sep*j)*normal_to_vrel[1]);
        fragment.vz = com.vz + frag_velocity*(cos(theta_sep*j)*unit_dvz + sin(theta_sep*j)*normal_to_vrel[2]);

        //Fragment radius is derived based on target's density
        fragment.r = get_radii(m_frag, rho);

        //Record collision
        fragment.last_collision = sim->t;


        //Add to mxsum vector to keep track of COM
        mxsum[0] +=fragment.m*fragment.x;
        mxsum[1] += fragment.m*fragment.y;    
        mxsum[2] += fragment.m*fragment.z;

        //Add to mvsum vector to keep track of momentum
        mvsum[0] += fragment.m*fragment.vx;
        mvsum[1] += fragment.m*fragment.vy;    
        mvsum[2] += fragment.m*fragment.vz;

        //Finally add fragment to simulation.
        reb_simulation_add(sim, fragment); 
    }

    //Now we correct for the COM and momentum offsets.
    //First, we need to see how much we are off from the COM in the begining. 
    double xoff[3] = {com.x - mxsum[0]/initial_mass, com.y - mxsum[1]/initial_mass, com.z - mxsum[2]/initial_mass};
    //Same for momentum
    double voff[3] = {com.vx - mvsum[0]/initial_mass, com.vy - mvsum[1]/initial_mass, com.vz - mvsum[2]/initial_mass};

    //Reassign new position and velocity to target (who is replaced by the largest remnant) to account for offsets:
    target -> x +=  xoff[0]*target->m/initial_mass;
    target -> y += xoff[1]*target->m/initial_mass; 
    target -> z += xoff[2]*target->m/initial_mass; 
    target -> vx += voff[0]*target->m/initial_mass; 
    target -> vy += voff[1]*target->m/initial_mass; 
    target -> vz += voff[2]*target->m/initial_mass; 

    //Reassign position and velocity to fragments to correct for offsets.
    for (int i = sim->N - n_frag; i < sim->N; i++){ 
        //mass fraction of fragment versus total initial mass
        double mass_fraction = sim->particles[i].m/initial_mass;

        sim->particles[i].x += xoff[0]*mass_fraction;
        sim->particles[i].y += xoff[1]*mass_fraction;
        sim->particles[i].z += xoff[2]*mass_fraction;

        sim->particles[i].vx += voff[0]*mass_fraction;
        sim->particles[i].vy += voff[1]*mass_fraction;
        sim->particles[i].vz += voff[2]*mass_fraction;
    }

    return 1; // Remove 1 particle from simulation (projectile)
}

