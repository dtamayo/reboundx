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

struct reb_simulation;
struct rebx_extras;
struct reb_particle;
struct reb_orbit;
struct reb_vec3d;
struct rebx_force;
struct rebx_operator;
enum REBX_COORDINATES;

void rebx_com_force(struct reb_simulation* const sim, struct rebx_force* const force, const enum REBX_COORDINATES coordinates, const int back_reactions_inclusive, const char* reference_name, struct reb_vec3d (*calculate_force) (struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p, struct reb_particle* source), struct reb_particle* const particles, const int N);

void rebx_tools_com_ptm(struct reb_simulation* const sim, struct rebx_operator* const operator, const enum REBX_COORDINATES coordinates, const int back_reactions_inclusive, const char* reference_name, struct reb_particle (*calculate_step) (struct reb_simulation* const sim, struct rebx_operator* const operator, struct reb_particle* p, struct reb_particle* source, const double dt), const double dt);

double rebx_Edot(struct reb_particle* const ps, const int N);

void rebx_calculate_jacobi_masses(const struct reb_particle* const ps, double* const m_j, const int N);

/****************************************
Effect helper functions
****************************************/
const double rebx_calculate_planet_trap(const double r, const double dedge, const double hedge);

// TLu 11/8/22
struct reb_vec3d rebx_tools_spin_and_orbital_angular_momentum(const struct rebx_extras* const rebx);
/*
void compute_transformation_angles(struct reb_simulation* sim, struct rebx_extras* rebx, double* theta1, double* theta2);
struct reb_vec3d EulerAnglesTransform(struct reb_vec3d xyz, const double Omega, const double I, const double omega);
void align_simulation(struct reb_simulation* sim, struct rebx_extras* rebx);
struct reb_vec3d rebx_transform_inv_to_planet(double inc, double omega, struct reb_vec3d spin_inv);
*/
#endif
