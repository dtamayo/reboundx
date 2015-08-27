#ifndef _LIBrebx_H
#define _LIBrebx_H
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

/*enum REBx_MODS {
	REBx_MODIFY_modify_orbits_forces 	= 0,
	REBx_MODIFY_modify_orbits_direct 	= 1,
	REBx_GR						= 2,
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

void rebx_add_modify_orbits_forces(struct reb_simulation* sim);
void rebx_add_modify_orbits_direct(struct reb_simulation* sim);
void rebx_add_gr(struct reb_simulation* sim, double c);

double* rebx_get_tau_a(struct reb_simulation* sim);
void rebx_set_tau_a(struct reb_simulation* sim, double* tau_a);

double* rebx_get_tau_e(struct reb_simulation* sim);
void rebx_set_tau_e(struct reb_simulation* sim, double* tau_e);

double* rebx_get_tau_inc(struct reb_simulation* sim);
void rebx_set_tau_inc(struct reb_simulation* sim, double* tau_inc);

double* rebx_get_tau_omega(struct reb_simulation* sim);
void rebx_set_tau_omega(struct reb_simulation* sim, double* tau_omega);

//void rebx_add(struct reb_simulation* sim, enum REBx_MODS perturbation);

/**
 * @cond PRIVATE
 * Internal functions used by reboundx.  User would not call these.
 */
void rebx_forces(struct reb_simulation* sim);
void rebx_ptm(struct reb_simulation* sim);

#endif
