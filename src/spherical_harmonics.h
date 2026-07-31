#ifndef spherical_harmonics_H
#define spherical_harmonics_H

#include "legendre.h"


/*
 * WARNING: rebx_spherical_harmonics_model is NOT thread-safe.
 *
 * The internal scratch buffers (P, dP, cosml, sinml) are owned by the model
 * and are overwritten on every call to rebx_spherical_harmonics_acceleration() /
 * rebx_spherical_harmonics_potential(). If this force loop is ever parallelized
 * (e.g. with OpenMP over particles or particle pairs), two threads calling
 * these functions on the SAME model concurrently will race on these buffers
 * and produce incorrect, non-deterministic results.
 *
 * Safe as-is for the current single-threaded, sequential-over-particles loop.
 * Do not parallelize any loop that calls these functions without first
 * giving each thread its own scratch buffers (or its own model instance).
 */


typedef struct {
    int n, m;
    double C, S;
} rebx_active_coeff;

typedef struct {        
    int N;

    double *A, *B;  
    double *D, *E;  
    double *F;      
    
    double *P, *dP;
    double *cosml, *sinml;
    double GM, R_eq;

    rebx_active_coeff *active;
    int n_active;
} rebx_spherical_harmonics_model;


/**
 * @brief Create a spherical harmonics (spherical_harmonics) gravity model.
 *
 * C and S are dense, 4pi fully-normalized coefficient arrays indexed by
 * idx(n,m) = n*(n+1)/2 + m, for 0 <= m <= n <= N (length (N+1)(N+2)/2
 * each). This dense, zero-padded layout is scanned once here to build an
 * internal sparse list of only the non-zero (n,m,C,S) terms; every
 * subsequent call to rebx_spherical_harmonics_acceleration() or
 * rebx_spherical_harmonics_potential() iterates over that sparse list instead of
 * the full dense triangle, so passing a large N with mostly-zero
 * coefficients costs almost nothing per timestep (the O(N^2) scan happens
 * only once, here). The effective working degree model->N is the highest
 * degree with a non-zero coefficient, and may be smaller than the input N.
 *
 * @param C     4pi fully-normalized cosine coefficients C_nm, dense
 *              triangular array indexed by idx(n,m). Dimensionless.
 * @param S     4pi fully-normalized sine coefficients S_nm, same layout
 *              and units as C.
 * @param GM    Gravitational parameter (G * mass) of the central body,
 *              in simulation units.
 * @param R_eq  Equatorial reference radius of the central body, in
 *              simulation length units.
 *
 * @return Newly allocated model, or NULL on allocation failure. Caller
 *         owns the returned pointer and must release it with
 *         rebx_spherical_harmonics_free().
 */



rebx_spherical_harmonics_model* rebx_spherical_harmonics_create(const double *C, const double *S, double GM, double R_eq);

/**
 * @brief Free all memory owned by a spherical_harmonics model. Safe on NULL.
 */

void rebx_spherical_harmonics_free(rebx_spherical_harmonics_model *model);


/**
 * @brief Evaluate the spherical_harmonics acceleration at a body-fixed position.
 *
 * NOT thread-safe: mutates the model's internal scratch buffers (P, dP,
 * cosml, sinml, rho) in place.
 *
 * @param model       spherical_harmonics model.
 * @param phi         Body-fixed latitude, radians.
 * @param r           Radial distance from the body's center, same length
 *                     units as R_eq.
 * @param lambda_body Body-fixed longitude, radians.
 * @param a_r         [out] Radial component of the acceleration.
 * @param a_phi       [out] Latitudinal component of the acceleration.
 * @param a_lambda    [out] Longitudinal component of the acceleration.
 * @param sinphi_out  [out, optional] sin(phi); NULL if not needed.
 * @param cosphi_out  [out, optional] cos(phi); NULL if not needed.
*/

void rebx_spherical_harmonics_acceleration(rebx_spherical_harmonics_model *model,
                                     double phi, double r,
                                     double lambda_body, 
                                     double *a_r, double *a_phi, double *a_lambda, double *sinphi_out, double *cosphi_out);


/**
 * @brief Evaluate the spherical_harmonics (gravitational potential) at a
 * body-fixed position.
 *
 * @param model       spherical_harmonics model.
 * @param phi         Body-fixed latitude, radians.
 * @param r           Radial distance from the body's center, same length
 *                     units as R_eq.
 * @param lambda_body Body-fixed longitude, radians.
 * @return Gravitational potential at (phi, r, lambda_body).
 */


double rebx_spherical_harmonics_potential(rebx_spherical_harmonics_model *model,
                                    double phi, double r, double lambda_body);
                            

#endif 
