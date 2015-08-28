/**
 * @file 	rebxtools.c
 * @brief 	Helper functions for reboundx
 * @author 	Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section 	LICENSE
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

#include <math.h>
#include <stdlib.h>
#include "rebound.h"
#include "rebxtools.h"

void rebxtools_orbit2p(double G, struct reb_particle* p, struct reb_particle* primary, struct reb_orbit o){
	int* err = malloc(sizeof(int));
	struct reb_particle p2 = reb_tools_orbit_to_particle_err(G,*primary, p->m, o.a, o.e, o.inc, o.Omega, o.omega, o.f, err);
	p->x = p2.x;
	p->y = p2.y;
	p->z = p2.z;
	p->vx = p2.vx;
	p->vy = p2.vy;
	p->vz = p2.vz;
}
