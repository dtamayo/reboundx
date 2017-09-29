/**
 * @file    integrator_implicit_midpoint.h
 * @brief   Interface for symplectic numerical integrator
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>, Hanno Rein
 *
 * @section     LICENSE
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
#ifndef _REBX_INTEGRATOR_IMPLICIT_MIDPOINT_H
#define _REBX_INTEGRATOR_IMPLICIT_MIDPOINT_H

void rebx_integrator_implicit_midpoint_integrate(struct reb_simulation* const sim, const double dt, struct rebx_effect* const effect);
#endif
