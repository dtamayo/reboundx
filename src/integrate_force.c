/**
 * @file    integrate_force.c
 * @brief   Generic operator for integrating a force across a timestep
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>, Hanno Rein
 *
 * @section LICENSE
 * Copyright (c) 2017 Dan Tamayo, Hanno Rein
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
#include <string.h>
#include <stdio.h>
#include "rebound.h"
#include "reboundx.h"
#include "core.h"
#include "rebxtools.h"

void rebx_integrate_force(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    struct rebx_extras* rebx = sim->extras;
    struct rebx_force* force = rebx_get_param(rebx, operator->ap, "force");
    if (force == NULL){
        reb_error(sim, "REBOUNDx Error: Force parameter not set in rebx_integrate operator. See examples for how to add as a parameter.\n");
    }
    enum rebx_integrator integrator = REBX_INTEGRATOR_EULER; // default
    enum rebx_integrator* integratorparam = rebx_get_param(rebx, operator->ap, "integrator");
    if (integratorparam != NULL){
        integrator = *integratorparam;
    }
    
    rebx_reset_accelerations(sim->particles, sim->N);

    switch(integrator){
        case REBX_INTEGRATOR_IMPLICIT_MIDPOINT:
        {
            rebx_integrator_implicit_midpoint_integrate(sim, dt, force);
            break;
        }
        case REBX_INTEGRATOR_RK2:
        {
            rebx_integrator_rk2_integrate(sim, dt, force);
            break;
        }
        case REBX_INTEGRATOR_RK4:
        {
            rebx_integrator_rk4_integrate(sim, dt, force);
            break;
        }
        case REBX_INTEGRATOR_EULER:
        {
            rebx_integrator_euler_integrate(sim, dt, force);
            break;
        }
        case REBX_INTEGRATOR_NONE:
        {
            break;
        }
        default:
        {
            break;
        }
    }
}
