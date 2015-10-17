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

enum REBX_PARAMS{
	ORB_TAU,
	TEST_INT,
	TEST_INT2,
};

enum REBX_COORDINATES{
	JACOBI,				// default.  REBX_COORDINATE members in structs get calloced to 0 and first enum always = 0
	BARYCENTRIC,
	HELIOCENTRIC
};

struct rebx_param{
	void* valPtr;
	enum REBX_PARAMS param;
	struct rebx_param* next;
	int rebx_index;
};

struct rebx_orb_tau{
	double tau_a;
	double tau_e;
	double tau_inc;
	double tau_omega;
	double tau_Omega;
};

struct rebx_param_to_be_freed{
	struct rebx_param_to_be_freed* next;
	struct rebx_param* param;
};

struct rebx_params_modify_orbits{
	double e_damping_p; // p paramseter from Deck & Batygin (2015) for how e-damping
	// is coupled to a-damping at order e^2
	// p = 1 : e-damping at const angular momentum.  p = 0 : no contribution to a-damping
	// equal to p/3 with p defined as in Goldreich & Schlichting 2014
	enum REBX_COORDINATES coordinates;
};

struct rebx_params_gr {
	double c;
};

struct rebx_extras {	
	struct reb_simulation* sim;
	struct rebx_param_to_be_freed** params_to_be_freed; // pointer to a linked list holding pointers to all
											// the allocated params for later freeing

	// these are pointers to simplify syntax.  Some structs need to update member variables
	// inside functions so we need to pass the pointer to them anyway

	struct rebx_params_modify_orbits* modify_orbits_forces;
	struct rebx_params_modify_orbits* modify_orbits_direct;
	struct rebx_params_gr* gr;

	// these are pointers to function pointers to use as arrays of function pointers for the user-added effects
	void (**ptm) (struct reb_simulation* const sim);
	void (**forces) (struct reb_simulation* const sim);

	int Nptm;
	int Nforces;
};

void* rebx_search_param(struct rebx_param* current, enum REBX_PARAMS param);
void rebx_add_param_orb_tau(struct reb_simulation* sim, struct rebx_param** paramsRef);
void rebx_set_tau_a(struct reb_simulation* sim, int p_index, double value);
double rebx_get_tau_a(struct reb_particle p);
void rebx_add_param_to_be_freed(struct rebx_extras* rebx, struct rebx_param* param);
void rebx_free_params(struct rebx_extras* rebx);
void rebx_free(struct rebx_extras* rebx);
void rebx_forces(struct reb_simulation* sim);
void rebx_ptm(struct reb_simulation* sim);
/*void rebx_set_double(struct reb_simulation* sim, int p_index, enum REBX_EXTRAS param, double value);
double rebx_get_double(struct reb_particle p, enum REBX_EXTRAS param);

void rebx_add_particle(struct rebx_extras* rebx, struct rebx_p_param* p_param);
void rebx_update_particles(struct rebx_extras* rebx, struct rebx_p_param* p_param);
void rebx_free_p_params(struct rebx_p_param* apPtr);
void rebx_free_particles(struct rebx_extras* rebx);
void rebx_free(struct rebx_extras* rebx);

void rebx_add_double_param(struct reb_simulation* sim, void** _p_paramsRef, enum REBX_EXTRAS param, double value);
void rebx_add_int_param(struct rebx_p_param** p_paramsRef, enum REBX_USER_PARAMS param, int value);

//double rebx_get_double_param(struct rebx_p_param*r p_param);
int rebx_get_int_param(struct rebx_p_param* p_param);
*/
struct rebx_extras* rebx_init(struct reb_simulation* sim);
void rebx_initialize(struct reb_simulation* sim, struct rebx_extras* rebx);

void rebx_add_modify_orbits_forces(struct reb_simulation* sim);
void rebx_add_modify_orbits_direct(struct reb_simulation* sim);
void rebx_add_gr(struct reb_simulation* sim, double c);
void rebx_add_gr_potential(struct reb_simulation* sim, double c);
void rebx_add_gr_implicit(struct reb_simulation* sim, double c);

#endif
