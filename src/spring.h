/**
 * @file        spring.h
 * @brief       Particle pair structure 
 * @author      Alice Quilen 
 */


#ifndef _SPRING_H
#define _SPRING_H


#include "rebound.h"

// spring structure 
struct rebx_spring  {
	double ks;   // spring constant
	double rs0;  // distance of no force
	double gamma; // damping coefficient
	double smax;  // 
	int    i;    // vertex 1 referring to a particle
	int    j;    // vertex 2 referring to a particle
};


extern double b_distance; // mush formation
extern double mush_distance; // mush spring connection 
extern double t_reform; // formation springs timescale
extern double gamma_all; // for gamma  of all springs
extern double Q_FORCERED; // force redistribution parameter

void spring_forces();
void zero_accel();

double spring_length();
double strain();
void springs_add();
void normalize ();
double mindist();
void centerbody();
void connect_springs_dist();
void rand_football();
void rand_football_from_sphere();
double Young_mush();
void set_gamma();
void spin();
void make_binary_spring();
void mom_inertia();
void measure_L();
void compute_semi();
void compute_semi_bin();
void total_mom();
double mean_L();
void spring_init(struct reb_simulation* r);
void output_png();
void output_png_single();
void spr_ang_mid();
void body_spin();
void print_tab();
void print_heat();
void invI();
double detI();
void eigenvalues();
void adjust_ks();
void adjust_mass_side();
void rotate_body();
double fill_hcp();
double fill_cubic();
double fill_cubic_cube();
double epsilon_mom();
double add_pluto_charon();
double add_one_mass();

#endif // _SPRING_H


