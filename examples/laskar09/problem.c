/**
 * Sun + Earth/Moon system example
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include <string.h>

double ss_pos[10][3] = 
{
	{0.0E-03  , 0.0E-04 , 0.0E-04},  
	{-1.36364695954795E-01 , 8.93397922857E-01 , 3.87458344639337E-01}, 

};
double ss_vel[10][3] = 
{
	{0.0E-06 ,  0.0E-06 ,  0.0E-08}, 
	{-1.73199988485296E-02,  -2.24430473176756E-03 ,  -9.73361115758044E-04 }, 
};

double ss_mass[10] =
{
	1.988544e30,
	1.988544e30/328900.53,
};

void heartbeat(struct reb_simulation* r);
double e_init;
double E0;
double tmax;
struct reb_vec3d_orb {
    float t; ///< time
    float energy; ///< semi-major axis
    float hamiltonian; ///< eccentricity
};

struct rebx_params_gr_potential* params;

//void lunar_quadrupole(struct reb_simulation* const r, struct rebx_effect* const effect);
void reb_output_binary_orbit(struct reb_simulation* r, FILE* of);
struct custom_params{
	double coefficient;
};

const double c = 173.26203208556151;
const int source_index = 0;

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();

	// Setup constants
	r->dt 			= 9.125;				// in days
	tmax			= 7.3e11;			// 2 Gyr
	r->G			= 1.4880826e-34;		// in AU^3 / kg / day^2.

	r->ri_whfast.safe_mode 	= 0;		// Turn off safe mode. Need to call reb_integrator_synchronize() before outputs. 
	r->ri_whfast.corrector 	= 11;		// 11th order symplectic corrector

	r->integrator		= REB_INTEGRATOR_WHFAST;
	r->heartbeat		= heartbeat;
	//r->integrator		= REB_INTEGRATOR_IAS15;		// Alternative non-symplectic integrator

	// Initial conditions
	for (int i=0;i<2;i++){
		struct reb_particle p = {0};
		p.x  = ss_pos[i][0]; 		p.y  = ss_pos[i][1];	 	p.z  = ss_pos[i][2];
		p.vx = ss_vel[i][0]; 		p.vy = ss_vel[i][1];	 	p.vz = ss_vel[i][2];
		p.m  = ss_mass[i];
		reb_add(r, p); 
	}

	reb_move_to_com(r);

	struct rebx_extras* rebx = rebx_init(r);  // first initialize rebx
// setting GR parameters
	rebx_add_gr_potential(rebx, source_index, c);

	params_gr = rebx_add_gr_potential(c=c, source_index=source_index)

// defining parameter for lunar quadrupole
//	const double G = r->G;
//	const double R_EM = 0.0025696*(1.0-(2.0/3.0)*(-3.9744e-13)*r->t);
//	const double em_factor = (-3.0/4.0)*G*R_EM*R_EM*0.9473;
//		
//	struct custom_params* params = malloc(sizeof(*params));
//	params->coefficient = em_factor;
//
//	int force_is_velocity_dependent = 0;
//	rebx_add_custom_force(rebx, lunar_quadrupole, force_is_velocity_dependent, params);

	e_init = reb_tools_energy(r);

	E0 = rebx_gr_potential_hamiltonian(r, params_gr); // relativistic hamiltonian

	system("rm -f relative_energy.txt");
	reb_integrate(r, tmax);
	rebx_free(rebx);                // Free all the memory allocated by rebx
}

void reb_output_binary_orbit(struct reb_simulation* r, FILE* of){
	const int N = r->N;
	if (of==NULL){
		reb_exit("Can not open file.");
	}
	struct reb_particle com = r->particles[0];
	for (int i=1;i<N;i++){
		struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com);
		struct reb_vec3d_orb porb;
		porb.t = r->t;
		porb.a = o.a;
		porb.e = o.e;
		fwrite(&(porb),sizeof(struct reb_vec3d_orb),1,of);
		com = reb_get_com_of_pair(com,r->particles[i]);
	}
}

	
void heartbeat(struct reb_simulation* r){

	if (reb_output_check(r,10*365.)){ // check if eccentricity is bigger than 0.4 every 10 years
		reb_output_timing(r, tmax);
		reb_integrator_synchronize(r);

		const int N = r->N;
		struct reb_particle com2 = r->particles[0];
		for (int i=1;i<N;i++){
			struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com2);
			if (o.e > 0.4 && o.e < 0.55){
				r->dt= 1.82625;
			}
			if (o.e > 0.55 && o.e < 0.70){
				r->dt= 1.0;
			}
			if (o.e > 0.70){
				r->dt = 1.0/4.0;
			}
			com2 = reb_get_com_of_pair(com2,r->particles[i]);

		}
	}

	if (reb_output_check(r, 1.e4*365.)){
		reb_integrator_synchronize(r);

		// relative energy
		FILE* f = fopen("/data_local/arachkov/laskar09/test/relative_energy.txt","a");
		double e = reb_tools_energy(r); 
		fprintf(f,"%e %e\t",r->t, fabs((e-e_init)/e_init));
		fclose(f);

		FILE* g = fopen("/data_local/arachkov/laskar09/test/hamiltonian_error3.txt","a");

		struct rebx_params_gr* params_gr = malloc(sizeof(*params_gr));
		params_gr->source_index = source_index;
		params_gr->c = c;

		double Ef = rebx_gr_potential_hamiltonian(r,params_gr); 
		fprintf(f,"%e\n",fabs((Ef-E0)/E0));
		fclose(g);

	}

}


void lunar_quadrupole(struct reb_simulation* const r, struct rebx_effect* const effect){ // lunar quadrupole correction on Earth, based on Quinn et al. 1991 paper
	
	// define constants from the r struct
	struct reb_particle* const particles = r->particles;

	const double dxsun = particles[3].x - particles[0].x;
	const double dysun = particles[3].y - particles[0].y;
	const double dzsun = particles[3].z - particles[0].z;

	const double r2sun = dxsun*dxsun + dysun*dysun + dzsun*dzsun;
	const double rsun = sqrt(r2sun);

	const double mass_sun = ss_mass[0];
	const double meml = 81.3007;
	const double massratio = 1.0/(meml + 2.0 + 1.0/meml);

	struct custom_params* params = effect->paramsPtr;
	double em_factor_earth = mass_sun * massratio * (params->coefficient)/(rsun*rsun*rsun*rsun*rsun);
	double em_factor_sun = 1./(mass_sun * massratio) * (params->coefficient)/(rsun*rsun*rsun*rsun*rsun);

	particles[3].ax += em_factor_earth*dxsun;
	particles[3].ay += em_factor_earth*dysun;
	particles[3].az += em_factor_earth*dzsun;

	particles[0].ax -= em_factor_sun*dxsun;
	particles[0].ay -= em_factor_sun*dysun;
	particles[0].az -= em_factor_sun*dzsun;
}
