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

void rebx_spring_forces(struct reb_simulation* const r, struct rebx_effect* const effect);
struct rebx_spring* rebx_get_param_springs(const void* const object, const char* const param_name);
struct rebx_spring* rebx_add_param_springs(void* object, uint32_t hash);
void rebx_set_param_spring(struct reb_simulation* const r, struct rebx_effect* effect, struct rebx_spring spr);
int add_spring_i(struct reb_simulation* const r, struct rebx_effect* effect, int i1, int i2, struct rebx_spring spring_vals);
void connect_springs_dist(struct reb_simulation* const r, struct rebx_effect* const effect, double h_dist, int i0, int imax, struct rebx_spring spring_vals);
double spring_length(struct reb_simulation* const r, struct rebx_spring spr);
int mym(int k);
void rand_football_from_sphere(struct reb_simulation* r, double dist, double ax, double by, double cz, double total_mass);
void rand_football(struct reb_simulation* const r, double dist, double ax, double by, double cz, double total_mass);
double fill_hcp(struct reb_simulation* r, double dd, double ax, double by, double cz, double total_mass);
double fill_cubic(struct reb_simulation* r, double dd, double ax, double by, double cz, double total_mass);
double fill_cubic_cube(struct reb_simulation* r, double dd, double ax, double by, double cz, double total_mass);
    
void zero_accel();

double strain();
void normalize ();
double mindist();
void centerbody();
void connect_springs_dist();
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
double epsilon_mom();
double add_pluto_charon();
double add_one_mass();

#endif // _SPRING_H


