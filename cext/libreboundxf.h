#ifndef _LIBREBOUNDXF_H
#define _LIBREBOUNDXF_H
#ifndef M_PI 
#define M_PI 3.1415926535879323846
#endif
/**
 * Particle structure.
 * @details This structure is used to represent one particle. Additional particle
 * properties should be added here. Note that when the data structure is changed, 
 * one must also update the equivalent declaration for MPI in communications_mpi.c.
 */
struct particle {
	double x;	/**< x-position of the particle. */
	double y;	/**< y-position of the particle. */
	double z;	/**< z-position of the particle. */
	double vx;	/**< x-velocity of the particle. */
	double vy;	/**< y-velocity of the particle. */
	double vz;	/**< z-velocity of the particle. */
	double ax;	/**< x-acceleration of the particle. */
	double ay;	/**< y-acceleration of the particle. */
	double az;	/**< z-acceleration of the particle. */
	double m;	/**< Mass of the particle. */
	int ID;		/**< Particle ID. */
};

/**
 * Struct representing a Keplerian orbit.
 */
struct orbit {
	double a;
	double r;	// Radial distance from central object
	double h;	// Angular momentum
	double P;	// Orbital period
	double l;
	double e;
	double inc;
	double Omega; 	// longitude of ascending node
	double omega; 	// argument of perihelion
	double f; 	// true anomaly
};

#endif
