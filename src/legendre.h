#ifndef LEGENDRE_H
#define LEGENDRE_H

static inline int idx(int n, int m)
{
    return n*(n+1)/2 + m;
}


/**
 * @brief Precompute recursion coefficients for the 4pi fully-normalized
 * associated Legendre function recursions. Call once per degree N before
 * compute_normalized_legendre() / compute_normalized_legendre_derivative().
 *
 * @param N Maximum degree.
 * @param A,B [out] size (N+1)(N+2)/2 — vertical recursion coefficients.
 * @param D,E [out] size N+1 — diagonal / sub-diagonal coefficients.
 * @param F   [out] size (N+1)(N+2)/2 — derivative recursion coefficient.
 */

void compute_recursion_coef(int N,
                            double *A, double *B,
                            double *D, double *E,
                            double *F);

/**
 * @brief Evaluate the 4pi fully-normalized associated Legendre functions
 * P_nm(sin(phi)) for all 0 <= m <= n <= N.
 *
 * @param xp sin(phi), phi being the body-fixed latitude, radians.
 * @param N  Maximum degree.
 * @param A,B,D,E Recursion coefficients from compute_recursion_coef().
 * @param P  [out] size (N+1)(N+2)/2, indexed by idx(n,m).
 */

void compute_normalized_legendre(double xp, int N,
                                 const double *A,
                                 const double *B,
                                 const double *D,
                                 const double *E,
                                 double *P);

/**
 * @brief Evaluate d/dphi of the normalized associated Legendre functions.
 * Requires P already filled by compute_normalized_legendre() for the same N.
 *
 * @param phi Body-fixed latitude, radians.
 * @param N   Maximum degree.
 * @param A,B,D,E,F Recursion coefficients from compute_recursion_coef().
 * @param P   Legendre function values, from compute_normalized_legendre().
 * @param dP  [out] size (N+1)(N+2)/2, indexed by idx(n,m).
 */

void compute_normalized_legendre_derivative(double phi, int N,
                                            const double *A,
                                            const double *B,
                                            const double *D,
                                            const double *E,
                                            const double *F,
                                            double *P,
                                            double *dP);

#endif

