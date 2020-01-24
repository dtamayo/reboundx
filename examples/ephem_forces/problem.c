/**
 * Highly eccentric orbits
 *
 * This example uses the IAS15 integrator to simulate
 * a very eccentric planetary orbit. The integrator
 * automatically adjusts the timestep so that the pericentre passages
 * resolved with high accuracy.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

double tmax;
void heartbeat(struct reb_simulation* r);
double jac0;

static void ephem(const int i, const double t, double* const m, double* const x, double* const y, double* const z){
    const double n = 1.;
    const double mu = 1.e-3;
    const double m0 = 1.-mu;
    const double m1 = mu; 
    if (i==0){
        *m = m0;
        const double mfac = -m1/(m0+m1);
        *x = mfac*cos(n*t);
        *y = mfac*sin(n*t);
        *z = 0.;
    }

    if (i==1){
        *m = m1;
        const double mfac = m0/(m0+m1);
        *x = mfac*cos(n*t);
        *y = mfac*sin(n*t);
        *z = 0.;
    }
}

double jac(struct reb_simulation* r){
    const double n = 1;
    const struct reb_particle* p = &r->particles[0];
    const double v2 = p->vx*p->vx + p->vy*p->vy;
    const double coriolis = 2*n*(p->x*p->vy-p->y*p->vx);
    double m, x, y, z;
    ephem(0, r->t, &m, &x, &y, &z);
    const double r1 = sqrt((x-p->x)*(x-p->x) + (y-p->y)*(y-p->y));
    const double mu1 = r->G*m;
    ephem(1, r->t, &m, &x, &y, &z);
    const double r2 = sqrt((x-p->x)*(x-p->x) + (y-p->y)*(y-p->y));
    const double mu2 = r->G*m;
    const double j = 2.*(mu1/r1 + mu2/r2) + coriolis - v2;
    return j;
}


int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->G            = 1;        // Gravitational constant
    r->integrator        = REB_INTEGRATOR_IAS15;
    r->heartbeat        = heartbeat;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->gravity = REB_GRAVITY_NONE;
    r->usleep = 10.;

    struct reb_particle tp = {0};
    tp.x = 1.5;
    tp.vy = 0.4;//1/sqrt(1.5); // for circular orbit

    reb_add(r, tp); 

    struct rebx_extras* rebx = rebx_attach(r);
    // Could also add "gr" or "gr_full" here.  See documentation for details.
    struct rebx_force* ephem_forces = rebx_load_force(rebx, "ephemeris_forces");
    rebx_add_force(rebx, ephem_forces); 
    // Have to set speed of light in right units (set by G & initial conditions).  Here we use default units of AU/(yr/2pi)
    rebx_set_param_int(rebx, &ephem_forces->ap, "N_ephem", 2);

    tmax            = 10000.;

    jac0 = jac(r);
    reb_integrate(r, tmax);
}

void heartbeat(struct reb_simulation* r){
    if(reb_output_check(r,tmax/10000.)){        // outputs to the screen
        reb_output_timing(r, tmax);
    }
    const double j = jac(r);
    fprintf(stderr, "jac const err: %e\n", fabs((j-jac0)/jac0));
}

