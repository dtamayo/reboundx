/**
 * @file    modify_orbits_forces.h
 * @brief   Update orbital elements with prescribed timescales using forces.
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
#ifndef REBX_MODIFY_ORBITS_FORCES_H
#define REBX_MODIFY_ORBITS_FORCES_H

struct reb_simulation;
struct rebx_effect;

void rebx_modify_orbits_forces(struct reb_simulation* const sim, struct rebx_effect* const effect);

#endif
