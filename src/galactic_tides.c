#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

/* Original defaults.
 *
 * Unit convention retained from your code:
 *   positions: AU
 *   masses: solar masses
 *   gt_rho0: solar masses / pc^3
 *   actual angular frequency: gt_omega0 * sqrt(sim->G)
 *
 * Changing length or mass units requires updating the conversion scheme.
 */

// Define 1pc in AU, according to IAU 2015 Resolution B2.
static const double GT_PC = 648000/M_PI;

// Define the default galactic coordinate axes
static const struct reb_vec3d GT_EX = {
    -0.0549427616647781,
    -0.993817643588233,
    -0.096476858538093
};
static const struct reb_vec3d GT_EY = {
     0.4941020309282433,
    -0.1110243460455739,
     0.8622857865103298
};
static const struct reb_vec3d GT_SGP = {
    8.6766610857711490e-01,
    2.9304928601984120e-04,
   -4.9714730025291526e-01
};

static double gt_get_double(
    struct rebx_extras* rebx,
    struct rebx_force* force,
    const char* name,
    const double fallback
){
    const double* value = rebx_get_param(rebx, force->ap, name);
    return value != NULL ? *value : fallback;
}

static struct reb_vec3d gt_get_vector(
    struct rebx_extras* rebx,
    struct rebx_force* force,
    const char* name,
    const struct reb_vec3d fallback
){
    const struct reb_vec3d* value =
        rebx_get_param(rebx, force->ap, name);

    return value != NULL ? *value : fallback;
}

void rebx_galactic_tides(
    struct reb_simulation* const sim,
    struct rebx_force* const tides,
    struct reb_particle* const particles,
    const int N
){
    struct rebx_extras* const rebx = sim->extras;

    // If parameters are not given, use the default values
    const double OMEGA =
        gt_get_double(rebx, tides, "gt_omega", 4.23201211e-9);
    const double RHO =
        gt_get_double(rebx, tides, "gt_rho", 0.1);
    const struct reb_vec3d EX =
        gt_get_vector(rebx, tides, "gt_ex", GT_EX);
    const struct reb_vec3d EY =
        gt_get_vector(rebx, tides, "gt_ey", GT_EY);
    const struct reb_vec3d SGP =
        gt_get_vector(rebx, tides, "gt_sgp", GT_SGP);

    // Determine if horizontal tides should be applied
    const int* horizontal_param =
        rebx_get_param(rebx, tides->ap, "gt_horizontal");

    const int horizontal =
        horizontal_param != NULL ? (*horizontal_param != 0) : 1;

    // Calculate the constant for the vertical tidal force
    const double vertical_coefficient =
        -4.0 * M_PI * sim->G * RHO / (GT_PC * GT_PC * GT_PC);

    double omega2 = 0.0;
    
    const int _N_real = sim->N - sim->N_var;

    if (horizontal){
        const double omega = OMEGA * sqrt(sim->G);
        const double phase = omega * sim->t;
        const double c = cos(phase);
        const double s = sin(phase);

        omega2 = omega * omega;

        // Rotate the horizontal basis vectors according to the current phase
        struct reb_vec3d ex = (struct reb_vec3d){EX.x * c + EY.x * s, EX.y * c + EY.y * s, EX.z * c + EY.z * s};
        struct reb_vec3d ey = (struct reb_vec3d){-EX.x * s + EY.x * c, -EX.y * s + EY.y * c, -EX.z * s + EY.z * c};
    }

    for (int i = sim->N_active; i < _N_real; i++){
        struct reb_particle* const p = &particles[i];

        // Calculate the galactic coordinates of the particle
        const double z_gal =
            SGP.x * p->x + SGP.y * p->y + SGP.z * p->z;

        // Calculate the vertical tidal forces
        const double tides_z = vertical_coefficient * z_gal;

        // Update particle accelerations for the vertical tide. 
        p->ax += SGP.x * tides_z;
        p->ay += SGP.y * tides_z;
        p->az += SGP.z * tides_z;

        if (horizontal){
            // Calculate the galactic coordinates of the particle in the horizontal plane
            const double x_gal = ex.x * p->x + ex.y * p->y + ex.z * p->z;
            const double y_gal = ey.x * p->x + ey.y * p->y + ey.z * p->z;

            // Calculate the horizontal tidal forces
            const double tides_x = omega2 * x_gal;
            const double tides_y = omega2 * y_gal;

            // Update particle accelerations if horizontal tides are considered
            p->ax += ex.x * tides_x + ey.x * tides_y;
            p->ay += ex.y * tides_x + ey.y * tides_y;
            p->az += ex.z * tides_x + ey.z * tides_y;
        }
    }
}