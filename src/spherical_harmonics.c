#include "spherical_harmonics.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

rebx_spherical_harmonics_model* rebx_spherical_harmonics_create(const double *C, const double *S, double GM, double R_eq) {

    rebx_spherical_harmonics_model *model = malloc(sizeof(*model));
    if (!model) return NULL;

    int count = 0;
    int Nreal = 1;
    int N = 50;
    for (int n = 2; n <= N; n++){
        for (int m = 0; m <= n; m++){
            const int k = idx(n,m);
            if (C[k] != 0.0 || S[k] != 0.0){
                count++;
                Nreal = n;
            }
        }
    }

    model->active = malloc(count * sizeof(rebx_active_coeff));
    if (!model->active) { rebx_spherical_harmonics_free(model); return NULL; }
    model->n_active = count;

    int a = 0;
    for (int n = 2; n <= N; n++){
        for (int m = 0; m <= n; m++) {
            const int k = idx(n, m);
            if (C[k] != 0.0 || S[k] != 0.0) {
                model->active[a].n = n;
                model->active[a].m = m;
                model->active[a].C = C[k];
                model->active[a].S = S[k];
                a++;
            }
        }
    }
    N = Nreal;
    model->N = N;
    const int size = (N + 1) * (N + 2) / 2;

    
    model->GM = GM;
    model->R_eq = R_eq;

    model->A     = calloc(size, sizeof(double));
    model->B     = calloc(size, sizeof(double));
    model->D     = calloc(N + 1, sizeof(double));
    model->E     = calloc(N + 1, sizeof(double));
    model->F     = calloc(size, sizeof(double));
    model->P     = calloc(size, sizeof(double));
    model->dP    = calloc(size, sizeof(double));
    model->cosml = calloc(N + 1, sizeof(double));
    model->sinml = calloc(N + 1, sizeof(double));

    if (!model->A || !model->B || !model->D || !model->E || !model->F ||
        !model->P || !model->dP || !model->cosml || !model->sinml) {
        rebx_spherical_harmonics_free(model);
        return NULL;
    }

    compute_recursion_coef(model->N, model->A, model->B, model->D, model->E, model->F);

    return model;
}

static inline void rebx_spherical_harmonics_compute_trig(rebx_spherical_harmonics_model *model, double lambda_body) {
    const double c = cos(lambda_body);
    const double s = sin(lambda_body);

    model->cosml[0] = 1.0;
    model->sinml[0] = 0.0;

    if (model->N >= 1) {
        model->cosml[1] = c;
        model->sinml[1] = s;
    }

    for (int m = 1; m < model->N; m++) {
        model->cosml[m+1] = c * model->cosml[m] - s * model->sinml[m];
        model->sinml[m+1] = s * model->cosml[m] + c * model->sinml[m];
    }
}


void rebx_spherical_harmonics_free(rebx_spherical_harmonics_model *model) {
    if (!model) return;
    free(model->A);
    free(model->B);
    free(model->D);
    free(model->E);
    free(model->F);
    free(model->P);
    free(model->dP);
    free(model->cosml);
    free(model->sinml);
    free(model->active);
    free(model);
}

/* NOT thread-safe: see warning in spherical_harmonics.h. Mutates model->P, model->dP,
 * model->cosml, model->sinml in place. */
void rebx_spherical_harmonics_acceleration(rebx_spherical_harmonics_model *model,
                                     double phi, double r,
                                     double lambda_body,
                                     double *a_r, double *a_phi, double *a_lambda, double *sinphi_out, double *cosphi_out) {

    const double sinphi = sin(phi);
    const double cosphi = cos(phi);
    const double GM = model->GM;
    const double R_eq = model->R_eq;

    compute_normalized_legendre(sinphi, model->N, model->A, model->B, model->D, model->E, model->P);
    compute_normalized_legendre_derivative(phi, model->N, model->A, model->B, model->D, model->E, model->F, model->P, model->dP);

    rebx_spherical_harmonics_compute_trig(model, lambda_body);

    double Vr = 0.0, Vphi = 0.0, Vlambda = 0.0;

    const double q = R_eq / r;
    double *rho = malloc((model->N + 1) * sizeof(double));
    rho[1] = q; 
    double acc = q * q;
    for (int n = 2; n <= model->N; n++) { rho[n] = acc; acc *= q; }

    for (int a = 0; a < model->n_active; a++) {
        const int n = model->active[a].n;
        const int m = model->active[a].m;
        const int k = idx(n, m);
        const double Pnm  = model->P[k];
        const double dPnm = model->dP[k];
        const double Cnm  = model->active[a].C;
        const double Snm  = model->active[a].S;

        const double trig  =  Cnm * model->cosml[m] + Snm * model->sinml[m];
        const double trig2 = -Cnm * model->sinml[m] + Snm * model->cosml[m];

        Vr      += (n + 1) * rho[n] * Pnm  * trig;
        Vphi    += rho[n] * dPnm * trig;
        Vlambda += rho[n] * m * Pnm * trig2;
    }

    free(rho);
    Vr      *= -GM / (r * r);
    Vphi    *=  GM / r;
    Vlambda *=  GM / r;

    *a_r      = Vr;
    *a_phi    = Vphi / r;
    *a_lambda = Vlambda / (r * cosphi);

    if (sinphi_out) *sinphi_out = sinphi;
    if (cosphi_out) *cosphi_out = cosphi;
}




static void rebx_spherical_harmonics_evaluate_basis(rebx_spherical_harmonics_model *model,
                                              double phi, double lambda_body,
                                              int need_derivative) {
    const double sinphi = sin(phi);
 
    compute_normalized_legendre(sinphi, model->N, model->A, model->B, model->D, model->E, model->P);
    if (need_derivative) {
        compute_normalized_legendre_derivative(phi, model->N, model->A, model->B, model->D, model->E, model->F, model->P, model->dP);
    }
    rebx_spherical_harmonics_compute_trig(model, lambda_body);
}

double rebx_spherical_harmonics_potential(rebx_spherical_harmonics_model *model,
                                    double phi, double r, double lambda_body) {
 
    const double GM   = model->GM;
    const double R_eq = model->R_eq;
 
    rebx_spherical_harmonics_evaluate_basis(model, phi, lambda_body, 0);
 
    double Vsum = 0.0;
    const double q = R_eq / r;
    double *rho = malloc((model->N + 1) * sizeof(double));
    rho[1] = q; 
    double acc = q * q;
    for (int n = 2; n <= model->N; n++) { rho[n] = acc; acc *= q; }

    for (int a = 0; a < model->n_active; a++) {
        const int n = model->active[a].n;
        const int m = model->active[a].m;
        const double Pnm  = model->P[idx(n, m)];
        const double Cnm  = model->active[a].C;
        const double Snm  = model->active[a].S;

        const double trig = Cnm * model->cosml[m] + Snm * model->sinml[m];
        Vsum += rho[n] * Pnm * trig;
    }
    free(rho);
    return (GM / r) * Vsum;
}


