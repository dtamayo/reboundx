/**
 * @file    rebxtools.h
 * @brief   Helper functions for reboundx
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

#ifndef _REBXTOOLS_H
#define _REBXTOOLS_H

#include "rebound.h"

#ifndef M_PI
#define M_PI 3.1415926535879323846 
#endif

struct reb_orbit rebxtools_particle_to_orbit_err(double G, struct reb_particle* p, struct reb_particle* primary, int* err);

struct reb_orbit rebxtools_particle_to_orbit(double G, struct reb_particle* p, struct reb_particle* primary);

void rebxtools_orbit2p(double G, struct reb_particle* p, struct reb_particle* primary, struct reb_orbit* o);

void rebxtools_orbit_to_particle(double G, struct reb_particle* p, struct reb_particle* primary, double a, double e, double inc, double Omega, double omega, double f, int* err);

void rebxtools_move_to_com(struct reb_simulation* const sim);

void rebxtools_update_com_with_particle(struct reb_particle* const com, const struct reb_particle* p);

void rebxtools_update_com_without_particle(struct reb_particle* const com, const struct reb_particle* p);

void rebxtools_get_com(const struct reb_simulation* const sim, const int first_N, struct reb_particle* com);

#endif
