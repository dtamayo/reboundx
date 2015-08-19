#ifndef _LIBrebxf_H
#define _LIBrebxf_H
#ifndef M_PI 
#define M_PI 3.1415926535879323846
#endif
#include "rebound.h"
#include "xftools.h"
#include "elements_direct.h"
#include "elements_forces.h"
#include "gr.h"

typedef void (*xfptr)(struct reb_simulation* const r);

/*enum REBXF_MODS {
	REBXF_MODIFY_ELEMENTS_FORCES 	= 0,
	REBXF_MODIFY_ELEMENTS_DIRECT 	= 1,
	REBXF_GR						= 2,
};*/

struct rebxf_mod_elements_forces {
	int allocatedN;
	double* tau_a;
	double* tau_e;
	double* tau_inc;
	double* tau_omega;
	double e_damping_p; // p parameter from Deck & Batygin (2015) for how e-damping
	// is coupled to a-damping at order e^2
	// p = 1 : e-damping at const angular momentum.  p = 0 : no contribution to a-damping
	// equal to p/3 with p defined as in Goldreich & Schlichting 2014
};

struct rebxf_mod_elements_direct {
	int allocatedN;
	double* tau_a;
	double* tau_e;
	double* tau_inc;
	double* tau_omega;
	double e_damping_p;
};

struct rebxf_mod_gr {
	int all_bodies;
};

struct rebxf_params {	
	struct reb_simulation* sim;
	xfptr* forces;
	xfptr* ptm;
	int Nforces;
	int Nptm;

	struct rebxf_mod_elements_forces elements_forces;
	struct rebxf_mod_elements_direct elements_direct;
	struct rebxf_mod_gr gr;
};

struct rebxf_params* rebxf_init(struct reb_simulation* sim);

void rebxf_add_elements_forces(struct reb_simulation* sim);
void rebxf_add_elements_direct(struct reb_simulation* sim);
void rebxf_add_gr(struct reb_simulation* sim);

double* rebxf_get_tau_a(struct reb_simulation* sim);
void rebxf_set_tau_a(struct reb_simulation* sim, double* tau_a);

double* rebxf_get_tau_e(struct reb_simulation* sim);
void rebxf_set_tau_e(struct reb_simulation* sim, double* tau_e);

double* rebxf_get_tau_inc(struct reb_simulation* sim);
void rebxf_set_tau_inc(struct reb_simulation* sim, double* tau_inc);

double* rebxf_get_tau_omega(struct reb_simulation* sim);
void rebxf_set_tau_omega(struct reb_simulation* sim, double* tau_omega);

struct rebxf_params* rebxf_init(struct reb_simulation* sim);
//void rebxf_add(struct reb_simulation* sim, enum REBXF_MODS perturbation);

/**
 * @cond PRIVATE
 * Internal functions used by reboundxf.  User would not call these.
 */
void rebxf_forces(struct reb_simulation* sim);
void rebxf_ptm(struct reb_simulation* sim);
void test(struct reb_simulation* sim);

#endif
