

// "fullo normalized associated Legendre functions" or 4 pi-normalized.




#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline int idx(int n, int m) {
    return n * (n + 1) / 2 + m;
}




void compute_recursion_coef(int N, double *A, double *B, double *D, double *E, double *F) {
    for (int n = 2; n <= N; n++) {
        for (int m = 0; m <= n - 2; m++) {

            double dn = n, dm = m; 

            A[idx(n, m)] = sqrt(((2 * dn + 1) * (2 * dn - 1)) / ((dn - dm) * (dn + dm)));

            B[idx(n, m)] = sqrt(((2 * dn + 1) * (dn + dm - 1) * (dn - dm - 1)) / ((2 * dn - 3) * (dn - dm) * (dn + dm)));
            
        }
    }

    for (int m = 1; m <= N; m++){
        D[m] = sqrt((2.0 * m + 1.0) / (2.0 * m));
    }

    for (int m = 0; m <= N - 1; m++) {
        E[m] = sqrt(2.0 * m + 3.0);
    }

    for (int m = 0; m <= N; m++) {
        for (int n = m + 1; n <= N; n++) {
            F[idx(n, m)] = sqrt((n - m) * (n + m + 1));
        }
    }
}



void compute_normalized_legendre(double xp, const int N, 
                                const double *A, const double *B, 
                                const double *D, const double *E, 
                                double *P) {

    const double u = sqrt(fmax(0.0, 1.0 - xp*xp));

    P[idx(0, 0)] = 1.0;

    // 1. Diagonal

    for (int m = 0; m <= N; m++) {
        if (m > 0)
            P[idx(m,m)] = -D[m] * u * P[idx(m-1,m-1)];
    // 2. Subdiagonal
        if (m < N)
            P[idx(m + 1, m)] = xp * E[m] * P[idx(m, m)];
    }

    // 3. Rest of elements
    for (int m = 0; m <= N - 2; m++) {
        for (int n = m + 2; n <= N; n++) {
            const int k = idx(n,m);
            P[k] = A[k] * xp * P[idx(n - 1, m)] - B[k] * P[idx(n - 2, m)];
        }
    }

}


void compute_normalized_legendre_derivative(double phi, const int N, 
                                const double *A, const double *B, 
                                const double *D, const double *E, 
                                const double *F,
                                double *P, double *dP) {

    const double xp = sin(phi);
    const double s = cos(phi);
    const double tan_phi = xp/s;

    for (int n = 0; n <= N; n++) {

        for (int m = 0; m < n; m++) {
            const int k = idx(n,m);
            dP[k] =
                m * tan_phi * P[k]
                - F[k] * P[idx(n,m+1)];

        }

        const int l = idx(n,n);
        dP[l] = n * tan_phi * P[l];
    }
}