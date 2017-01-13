/** * Solar System */ 
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include <string.h>
#include <time.h>

double ss_pos[10][3] = 
{
	{0.0E-03  , 0.0E-04 , 0.0E-04},  
//	{3.43926450169642E-01 +1*62.5*0.38*6.68459e-15, 4.56154799533135E-02 +1*62.5*0.38*6.68459e-15, -1.09240372119325E-02 +1*62.5*0.38*6.68459e-15}, 
//	{1.42965184343246E-01 , 6.47005066033887E-01  , 2.8248240066038E-01 }, 
	{-1.36364695954795E-01 , 8.93397922857E-01 , 3.87458344639337E-01}, 
//	{-1.36983397618342E-00  , 8.43135248017904E-01  , 4.2383290661143E-01 }, 
//	{3.34936422369601E+00 , -3.47376144901258E+00  , -1.5721496863938E-00 }, 
//	{-8.97250506828211E+00 , 2.27968200813286E-00 , 1.33033860971146E-00 }, 
//	{-1.00300399532732E+00  , 1.73235084732718E+01  , 7.60482504641591E-00},  
//	{-2.91945807386112E+01  , -7.71928519199847E+00 , -2.42724537828877E-00},  
//	{-2.62336507820155E+01  , 2.05619754200559E+01 , 1.44445571277807E+01 }, 

};

double ss_vel[10][3] = 
{
	{0.0E-06 ,  0.0E-06 ,  0.0E-08}, 
//	{-8.46653204500986E-03,  2.56146053072818E-02,  1.45868453362396E-02 }, 
//	{-1.98938122793425E-02,  3.11311946611859E-03,  2.6594458057458E-03 },
	{-1.73199988485296E-02,  -2.24430473176756E-03 ,  -9.73361115758044E-04 }, 
//	{-7.38456383127117E-03 ,  -9.4773586392742E-03,  -4.15165513666213E-03 }, 
//	{5.58564314958231E-03,  4.96226113722251E-03,  1.99227692673937E-03 }, 
//	{-1.85825195699671E-03 ,  -4.98385858141744E-03,  -1.98025741280725E-03}, 
//	{-3.95525416301772E-03,  -3.75913785112941E-04 ,  -1.08849991287794E-04 }, 
//	{8.20748057818085E-04 ,  -2.77209825958023E-03 ,  -1.15611603043592E-03}, 
//	{-1.31588869828116E-03 ,  -2.62012820549352E-03 ,  -4.27083355026298E-04},
};

double ss_mass[10] =
{
	1.988544e30,
//	1.988544e30/6023600,
//	1.988544e30/408523.5,
	1.988544e30/328900.53,
//	1.988544e30/3098710,
//	1.988544e30/1047.355,
//	1.988544e30/3498.5,
//	1.988544e30/22869,
//	1.988544e30/19314,
//	1.988544e30/3000000,
};

struct reb_vector {
    double t; ///< time
    double dt; ///< timestep
    double wall_time; // time passed to compute
    double energy; // tracking energy
    double angular_momentum; // tracking angular momentum
//    double linear_momentum; // tracking linear momentum
};

struct reb_cartesian_vector {
    double x; 
    double y; 
    double z; 
    double vx; 
    double vy; 
    double vz; 
};
/*
// lunar quadrupole force correction
void lunar_quadrupole(struct reb_simulation* const r, struct rebx_effect* const effect);
// calculate energy from lunar correction
double lunar_energy(struct reb_simulation* const r);

struct custom_params{
	double coefficient;
};
*/
double E0; // initial energy
//struct rebx_params_gr_potential* params_gr;
double L0; // initial angular momentum
clock_t begin; // beginning time
const double c = 173.26203208556151; //speed of light in AU/day
const int source_index = 0;
void heartbeat(struct reb_simulation* r);
double tmax;

// function to output data in binary
void reb_output_binary_orbit(struct reb_simulation* r, char* filename);

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();

	// Setup constants
	r->dt 			= 9.125;				// in days
	tmax			= 7.3e8;			
	//tmax			= 7.3e11;			// 2 Gyr
	r->G			= 1.4880826e-34;		// in AU^3 / kg / day^2.

	//r->ri_whfast.safe_mode 	= 0;		// Turn off safe mode. Need to call reb_integrator_synchronize() before outputs. 
	//r->ri_whfast.corrector 	= 11;		// 11th order symplectic corrector

	//r->integrator		= REB_INTEGRATOR_WHFAST;
	r->heartbeat		= heartbeat;
	r->integrator		= REB_INTEGRATOR_IAS15;		// Alternative non-symplectic integrator

	// Initial conditions
	for (int i=0;i<2;i++){
	//for (int i=0;i<10;i++){
		struct reb_particle p = {0};
		p.x  = ss_pos[i][0]; 		p.y  = ss_pos[i][1];	 	p.z  = ss_pos[i][2];
		p.vx = ss_vel[i][0]; 		p.vy = ss_vel[i][1];	 	p.vz = ss_vel[i][2];
		p.m  = ss_mass[i];
		reb_add(r, p); 
	}

	reb_move_to_com(r);
/*
	struct rebx_extras* rebx = rebx_init(r);  // first initialize rebx
        params_gr = rebx_add_gr_potential(rebx, source_index, c);
        E0 = rebx_gr_potential_hamiltonian(r, params_gr); // relativistic hamiltonian

// defining parameter for lunar quadrupole
	double G = r->G;
	double R_EM = 0.0025696*pow((1.0+(0.101773133860118/0.113381622646105)*1.e-9*(r->t/365.)-0.0178457555910623*1.e-18*(r->t/365.)*(r->t/365.)),0.113381622646105);
	double em_factor = (-3.0/4.0)*G*R_EM*R_EM*0.9473;
		
	struct custom_params* params = malloc(sizeof(*params));
	params->coefficient = em_factor;

	int force_is_velocity_dependent = 0;
	rebx_add_custom_force(rebx, lunar_quadrupole, force_is_velocity_dependent, params);
*/		
//	E0 = reb_tools_energy(r) + lunar_energy(r);
	E0 = reb_tools_energy(r);
	struct reb_vec3d L = reb_tools_angular_momentum(r);
	L0 = sqrt(L.x*L.x + L.y*L.y + L.z*L.z);
	system("rm -f test.bin");

	begin = clock();

	reb_integrate(r, tmax);
//	rebx_free(rebx);                // Free all the memory allocated by rebx
}

void heartbeat(struct reb_simulation* r){

	if (reb_output_check(r, 1.e0*365.)){
		reb_output_timing(r, tmax);
		//reb_integrator_synchronize(r);
		char* of = "test.bin";
		reb_output_binary_orbit(r,of);
	}

}

void reb_output_binary_orbit(struct reb_simulation* r, char* filename){
	FILE* of = fopen(filename,"ab"); 
	const int N = r->N;
	if (of==NULL){
		reb_exit("Can not open file.");
	}
//	struct reb_particle com = r->particles[0];

	struct reb_vec3d L = reb_tools_angular_momentum(r);
	double Lf = sqrt(L.x*L.x + L.y*L.y + L.z*L.z); // absolute value of total angular momentum

//	double Ef = reb_tools_energy(r) + lunar_energy(r); // calculating total energy
	double Ef = reb_tools_energy(r); // calculating total energy
//      double Ef = rebx_gr_potential_hamiltonian(r, params_gr); // relativistic hamiltonian

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC; // calculating wall time

	struct reb_vector vec;
	vec.t = r->t;
	vec.dt = r->dt;
	vec.wall_time = time_spent;
	vec.energy = (Ef-E0)/E0;
	vec.angular_momentum = (Lf-L0)/L0;
	fwrite(&(vec),sizeof(struct reb_vector),1,of);

	for (int i=1;i<N;i++){
//		struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com);

		struct reb_cartesian_vector vec_car;
		vec_car.x = r->particles[i].x; 
		vec_car.y = r->particles[i].y; 
		vec_car.z = r->particles[i].z; 
		vec_car.vx = r->particles[i].vx; 
		vec_car.vy = r->particles[i].vy; 
		vec_car.vz = r->particles[i].vz; 
		fwrite(&(vec_car),sizeof(struct reb_cartesian_vector),1,of);
	}
	fclose(of);
}
/*
void lunar_quadrupole(struct reb_simulation* const r, struct rebx_effect* const effect){ // lunar quadrupole correction on Earth, based on Quinn et al. 1991 paper
	
	// define constants from the r struct
	struct reb_particle* const particles = r->particles;

	double dxsun = particles[1].x - particles[0].x;
	double dysun = particles[1].y - particles[0].y;
	double dzsun = particles[1].z - particles[0].z;

	double r2sun = dxsun*dxsun + dysun*dysun + dzsun*dzsun;
	double rsun = sqrt(r2sun);

	double mass_sun = particles[0].m;
	double meml = 81.3007; // Earth-to-moon mass ratio
	double massratio = 1.0/(meml + 2.0 + 1.0/meml);

	struct custom_params* params = effect->paramsPtr;
	double em_factor_earth = mass_sun * massratio * (params->coefficient)/(rsun*rsun*rsun*rsun*rsun);
	
	particles[1].ax += em_factor_earth*dxsun;
	particles[1].ay += em_factor_earth*dysun;
	particles[1].az += em_factor_earth*dzsun;

	particles[0].ax -= (particles[1].m/particles[0].m)*em_factor_earth*dxsun;
	particles[0].ay -= (particles[1].m/particles[0].m)*em_factor_earth*dysun;
	particles[0].az -= (particles[1].m/particles[0].m)*em_factor_earth*dzsun;
}

double lunar_energy(struct reb_simulation* const r){
	
	// define constants from the r struct
	struct reb_particle* const particles = r->particles;

	double dxsun = particles[1].x - particles[0].x;
	double dysun = particles[1].y - particles[0].y;
	double dzsun = particles[1].z - particles[0].z;

	double r2sun = dxsun*dxsun + dysun*dysun + dzsun*dzsun;
	double rsun = sqrt(r2sun);

	double mass_sun = particles[0].m;
	double meml = 81.3007;
	double massratio = 1.0/(meml + 2.0 + 1.0/meml);

	double R_EM = 0.0025696*pow((1.0+(0.101773133860118/0.113381622646105)*1.e-9*(r->t/365.)-0.0178457555910623*1.e-18*(r->t/365.)*(r->t/365.)),0.113381622646105);
	double em_factor = (-3.0/4.0)*r->G*R_EM*R_EM*0.9473;
	double EM_energy = (1./3.)*r->particles[1].m * mass_sun * massratio * (em_factor)/(rsun*rsun*rsun);
	return EM_energy;
}*/
