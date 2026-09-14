/**
 * constants.h
 * Introduces the constants structure from Earth and other celestial bodies. This structure is used to pass the constants to the C code.
 * It's worth pointing out that the lists C and S are the unnormalized coefficients up to degree 20 from the EGM2008 model for the Earth 
 * and the AIUB-GRL350A model for the Moon.
 * This models are available at the ICGEM webside: https://icgem.gfz.de/tom_celestial.
 * While the other constants (in the International System of Units) are available at the IERS Conventions (2010) and the IAU WGCCRE report (Archinal et al. 2018). 
 * Or, for instance:  Vallado, Fundamentals of Astrodynamics and Applications.
 */ 



#ifndef REBX_CONSTANTS_H
#define REBX_CONSTANTS_H

typedef struct {
    double M;
    double R_eq;
    double omega;
    double theta0;
    double epoch;
    const double *C_coeffs;
    const double *S_coeffs;
} rebx_gravity_constants;

#endif