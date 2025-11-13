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
 * $Collisions$     // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script). 
 * 
 * ======================= ===============================================
 * Authors                 H. Rein 
 * Based on                None
 * C Example               :ref:`c_example_merging_collisions`
 * ======================= ===============================================
 * 
 * This is a simple example implementation of a REBOUNDx collision module.
 * The outcome is similar to the built-in REBOUND function reb_collision_resolve_merge.
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

enum REB_COLLISION_RESOLVE_OUTCOME rebx_merging_collisions(struct reb_simulation* const sim, struct rebx_collision_resolve* const collision_resolve, struct reb_collision c){
    struct reb_particle* pi = &(sim->particles[c.p1]);
    struct reb_particle* pj = &(sim->particles[c.p2]);

    double invmass = 1.0/(pi->m + pj->m);

    // Merge by conserving mass, volume and momentum
    pi->vx = (pi->vx*pi->m + pj->vx*pj->m)*invmass;
    pi->vy = (pi->vy*pi->m + pj->vy*pj->m)*invmass;
    pi->vz = (pi->vz*pi->m + pj->vz*pj->m)*invmass;
    pi->x  = (pi->x*pi->m + pj->x*pj->m)*invmass;
    pi->y  = (pi->y*pi->m + pj->y*pj->m)*invmass;
    pi->z  = (pi->z*pi->m + pj->z*pj->m)*invmass;
    pi->m  = pi->m + pj->m;
    pi->r  = cbrt(pi->r*pi->r*pi->r + pj->r*pj->r*pj->r);
    pi->last_collision = sim->t;

    return REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P2; // Remove 2 particle from simulation
}

