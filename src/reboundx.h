/**
 * @file 	reboundx.h
 * @brief 	REBOUNDx API definition.
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

#ifndef _LIBREBX_H
#define _LIBREBX_H
#ifndef M_PI 
#define M_PI 3.1415926535879323846
#endif
#define C_DEFAULT 10064.9150404 // speed of light in default units of AU/(yr/2pi)

#include "rebound.h"
#include "rebxtools.h"
#include "modify_orbits_direct.h"
#include "modify_orbits_forces.h"
#include "gr.h"

typedef void (*xptr)(struct reb_simulation* const r);

/*enum REBX_EXTRAS {
	REBX_MODIFY_ORBITS_FORCES	= 0,
	REBX_MODIFY_ORBITS_DIRECT	= 1,
	REBX_GR						= 2,
};*/

struct rebx_params_modify_orbits_forces {
	int allocatedN;
	double* tau_a;
	double* tau_e;
	double* tau_inc;
	double* tau_omega;
	double e_damping_p; // p paramseter from Deck & Batygin (2015) for how e-damping
	// is coupled to a-damping at order e^2
	// p = 1 : e-damping at const angular momentum.  p = 0 : no contribution to a-damping
	// equal to p/3 with p defined as in Goldreich & Schlichting 2014
};

struct rebx_params_modify_orbits_direct {
	int allocatedN;
	double* tau_a;
	double* tau_e;
	double* tau_inc;
	double* tau_omega;
	double e_damping_p;
};

struct rebx_params_gr {
	double c;
};

struct rebx_extras {	
	struct reb_simulation* sim;
	xptr* forces;
	xptr* ptm;
	int Nforces;
	int Nptm;

	struct rebx_params_modify_orbits_forces modify_orbits_forces;
	struct rebx_params_modify_orbits_direct modify_orbits_direct;
	struct rebx_params_gr gr;
};

struct rebx_extras* rebx_init(struct reb_simulation* sim);
void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx);

void rebx_add_modify_orbits_forces(struct reb_simulation* sim);
void rebx_add_modify_orbits_direct(struct reb_simulation* sim);
void rebx_add_gr(struct reb_simulation* sim, double c);
void rebx_add_gr_potential(struct reb_simulation* sim, double c);
void rebx_add_gr_implicit(struct reb_simulation* sim, double c);

double* rebx_get_tau_a(struct reb_simulation* sim);
void rebx_set_tau_a(struct reb_simulation* sim, double* tau_a);

double* rebx_get_tau_e(struct reb_simulation* sim);
void rebx_set_tau_e(struct reb_simulation* sim, double* tau_e);

double* rebx_get_tau_inc(struct reb_simulation* sim);
void rebx_set_tau_inc(struct reb_simulation* sim, double* tau_inc);

double* rebx_get_tau_omega(struct reb_simulation* sim);
void rebx_set_tau_omega(struct reb_simulation* sim, double* tau_omega);

//void rebx_add(struct reb_simulation* sim, enum REBX_EXTRAS extra);

/**
 * @cond PRIVATE
 * Internal functions used by reboundx.  User would not call these.
 */
void rebx_forces(struct reb_simulation* sim);
void rebx_ptm(struct reb_simulation* sim);

void rebx_free_extras(struct rebx_extras* const rebx);
void rebx_free_pointers(struct rebx_extras* const rebx);

#endif
