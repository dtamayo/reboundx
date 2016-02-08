/**
 * @file    reboundx.c
 * @brief   Initialization and getter/setter routines.
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

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include "core.h"
#include "reboundx.h"

const char* rebx_build_str = __DATE__ " " __TIME__; // Date and time build string. 
const char* rebx_version_str = "2.5.0";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.

/****************************************
User API
*****************************************/

struct rebx_extras* rebx_init(struct reb_simulation* sim){
    struct rebx_extras* rebx = malloc(sizeof(*rebx));
    rebx_initialize(sim, rebx);
    return rebx;
}

void rebx_free(struct rebx_extras* rebx){
    rebx_free_pointers(rebx);
    free(rebx);
}

void rebx_add_custom_post_timestep_modification(struct rebx_extras* rebx, void (*custom_post_timestep_modification)(struct reb_simulation* const sim, struct rebx_effect* custom_effect), void* custom_params){
    struct reb_simulation* sim = rebx->sim;
    sim->post_timestep_modifications = rebx_post_timestep_modifications;

    rebx_add_post_timestep_modification(rebx, custom_params, "custom_post_timestep_modification", custom_post_timestep_modification);
}

void rebx_add_custom_force(struct rebx_extras* rebx, void (*custom_force)(struct reb_simulation* const sim, struct rebx_effect* custom_effect), int force_is_velocity_dependent, void* custom_params){
    struct reb_simulation* sim = rebx->sim;
    sim->additional_forces = rebx_forces;
    if (force_is_velocity_dependent){
        sim->force_is_velocity_dependent = 1;
    }
    rebx_add_force(rebx, custom_params, "custom_force", custom_force);
}

/* Function to test whether REBOUNDx can load librebound.so and call REBOUND functions. */

double install_test(void){
    struct reb_simulation* sim = reb_create_simulation();
    struct reb_particle p = {0};
    p.m = 1.; 
    reb_add(sim, p); 
    struct reb_particle p1 = reb_tools_orbit2d_to_particle(sim->G, p, 0., 1., 0.2, 0., 0.);
    reb_add(sim, p1);
    reb_integrate(sim, 1.);
    return sim->particles[1].x;
}
