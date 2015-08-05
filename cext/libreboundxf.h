#ifndef _LIBREBOUNDXF_H
#define _LIBREBOUNDXF_H
#ifndef M_PI 
#define M_PI 3.1415926535879323846
#endif
#include "rebound.h"

void reboundxf_forces(struct reb_simulation* const sim);

void reboundxf_modify_elements(struct reb_simulation* const sim);

double reboundxf_get_e_damping_p(void);
void reboundxf_set_e_damping_p(double val);

void reboundxf_set_migration(double *_tau_a, int N);
void reboundxf_set_e_damping(double *_tau_e, int N);
void reboundxf_set_i_damping(double *_tau_i, int N);
void reboundxf_set_peri_precession(double *_tau_po, int N);

/**
 * Structure representing one REBOUND particle.
 */
/*struct reb_particle {
	double x;			///< x-position of the particle. 
	double y;			///< y-position of the particle. 
	double z;			///< z-position of the particle. 
	double vx;			///< x-velocity of the particle. 
	double vy;			///< y-velocity of the particle. 
	double vz;			///< z-velocity of the particle. 
	double ax;			///< x-acceleration of the particle. 
	double ay;			///< y-acceleration of the particle. 
	double az;			///< z-acceleration of the particle. 
	double m;			///< Mass of the particle. 
	double r; 			///< Radius of the particle. 
	double lastcollision;		///< Last time the particle had a physical collision.
	struct reb_treecell* c;		///< Pointer to the cell the particle is currently in.
	int id;				///< Unique id to identify particle.
};

/**
 * Struct representing a Keplerian orbit.
 */
/*struct reb_orbit {
	double a;	///< Semi-major axis
	double r;	///< Radial distance from central object
	double h;	///< Angular momentum
	double P;	///< Orbital period
	double l;	///< Mean longitude
	double e;	///< Eccentricity
	double inc;	///< Inclination
	double Omega; 	///< Longitude of ascending node
	double omega; 	///< Argument of perihelion
	double f; 	///< True anomaly
};*/

#endif
