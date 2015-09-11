/**
 * @file 	gr.h
 * @brief 	Post-newtonian general relativity corrections
 * @author 	Pengshuai (Sam) Shi <tamayo.daniel@gmail.com>
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Pengshuai (Sam) Shi, Dan Tamayo, Hanno Rein
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

#ifndef _GR_H
#define _GR_H

#include "rebound.h"

void rebx_gr(struct reb_simulation* const sim);
void rebx_gr_potential(struct reb_simulation* const sim);
void rebx_gr_single_mass(struct reb_simulation* const sim);

#endif
