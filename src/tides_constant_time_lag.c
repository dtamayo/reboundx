/**
 * @file    tides_constant_time_lag.c
 * @brief   Add constant time lag tides raised on primary, orbiting bodies, or both
 * @author  Stanley A. Baronett <stanley.a.baronett@gmail.com>, Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Tides$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 Stanley A. Baronett, D. Tamayo, Noah Ferich
 * Implementation Paper    Baronett et al., in prep.
 * Based on                `Hut 1981 <https://ui.adsabs.harvard.edu/#abs/1981A&A....99..126H/abstract>`_, `Bolmont et al., 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...583A.116B/abstract>`_.
 * C Example               :ref:`c_example_tides_constant_time_lag`.
 * Python Example          `TidesConstantTimeLag.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/TidesConstantTimeLag.ipynb>`_.
 * ======================= ===============================================
 *
 * This adds constant time lag tidal interactions between orbiting bodies in the simulation and the primary, both from tides raised on the primary and on the other bodies.
 * In all cases, we need to set masses for all the particles that will feel these tidal forces. After that, we can choose to include tides raised on the primary, on the "planets", or both, by setting the respective bodies' physical radius particles[i].r, k1 (apsidal motion constant, half the tidal Love number), constant time lag tau, and rotation rate Omega. See Hut (1981) and Bolmont et al. 2015 above.
 *
 * If tau is not set, it will default to zero and yield the conservative piece of the tidal potential.
 * 
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * particles[i].r (float)       Yes         Physical radius (required for contribution from tides raised on the body).
 * tctl_k1 (float)                   Yes         Apsidal motion constant (half the tidal Love number k2).
 * tctl_tau (float)                  No          Constant time lag. If not set will default to 0 and give conservative tidal potential
 * Omega (float)                No          Rotation rate. If not set will default to 0
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "reboundx.h"

static void rebx_calculate_tides(struct reb_particle* source, struct reb_particle* target, const double G, const double k1, const double tau, const double Omega){
    const double ms = source->m;
    const double mt = target->m;
    const double Rt = target->r;

    const double mratio = ms/mt; // have already checked for 0 and inf
    const double fac = mratio*k1*Rt*Rt*Rt*Rt*Rt; 
    
    const double dx = target->x - source->x; 
    const double dy = target->y - source->y;
    const double dz = target->z - source->z;
    const double dr2 = dx*dx + dy*dy + dz*dz; 
    const double prefac = -3*G/(dr2*dr2*dr2*dr2)*fac;
    double rfac = prefac;

    if (tau != 0){
        const double dvx = target->vx - source->vx; 
        const double dvy = target->vy - source->vy;
        const double dvz = target->vz - source->vz;

        rfac *= (1. + 3.*tau/dr2*(dx*dvx + dy*dvy + dz*dvz));
        const double thetafac = -prefac*tau;

        const double hx = dy*dvz - dz*dvy;
        const double hy = dz*dvx - dx*dvz;
        const double hz = dx*dvy - dy*dvx;

        const double thetadotcrossrx = (hy*dz - hz*dy)/dr2; // vec(thetadot) cross vec(r) Bolmont Eq. 7
        const double thetadotcrossry = (hz*dx - hx*dz)/dr2;
        const double thetadotcrossrz = (hx*dy - hy*dx)/dr2;

        // Assumes all spins vec(Omega) = Omega zhat, i.e., spin fixed along z axis
        const double Omegacrossrx = -Omega*dy;
        const double Omegacrossry = Omega*dx; 
        const double Omegacrossrz = 0.;

        target->ax += thetafac*ms*(Omegacrossrx-thetadotcrossrx);
        target->ay += thetafac*ms*(Omegacrossry-thetadotcrossry);
        target->az += thetafac*ms*(Omegacrossrz-thetadotcrossrz);
        source->ax -= thetafac*mt*(Omegacrossrx-thetadotcrossrx);
        source->ay -= thetafac*mt*(Omegacrossry-thetadotcrossry);
        source->az -= thetafac*mt*(Omegacrossrz-thetadotcrossrz);
    }

    target->ax += rfac*ms*dx;
    target->ay += rfac*ms*dy;
    target->az += rfac*ms*dz;
    source->ax -= rfac*mt*dx;
    source->ay -= rfac*mt*dy;
    source->az -= rfac*mt*dz;
}


void rebx_tides_constant_time_lag(struct reb_simulation* const sim, struct rebx_force* const tides, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    const double G = sim->G;

    struct reb_particle* star = &particles[0];  // assumes nearly Keplerian motion around a single primary (particles[0])
   
    double starmass = star->m;
    if (star->m == 0){                          // nothing makes sense if primary has no mass
        return;
    }
    
    double starradius = star->r;
    double starOmega = 0.;
    double stark1 = 0.; 
    double fixedstartau = 0.;
    double* starOmegaptr = rebx_get_param(rebx, star->ap, "Omega");
    double* stark1ptr = rebx_get_param(rebx, star->ap, "tctl_k1");
    double* fixedstartauptr = rebx_get_param(rebx, star->ap, "tctl_tau");
    if (stark1ptr){
        stark1 = *stark1ptr;
    }
    if (starOmegaptr){
        starOmega = *starOmegaptr;
    }
    if (fixedstartauptr){
        fixedstartau = *fixedstartauptr;
    }

    for (int i=1; i<N; i++){
        struct reb_particle* planet = &particles[i]; 
    
        // Calculate tides raised on the planets
        double* k1 = rebx_get_param(rebx, planet->ap, "tctl_k1");
        if (k1 != NULL && planet->r > 0 && planet->m > 0){ // tides only nonzero if k1 and finite size & mass are set
            double tau = 0.;
            double Omega = 0.;
            double* tauptr = rebx_get_param(rebx, planet->ap, "tctl_tau");
            // We don't require time lag tau to be set. Might just want conservative piece of tidal potential
            if (tauptr){
                tau = *tauptr;
                double* Omegaptr = rebx_get_param(rebx, planet->ap, "Omega");
                if (Omegaptr){
                    Omega = *Omegaptr;
                }
            }
            rebx_calculate_tides(star, planet, G, *k1, tau, Omega);
        }
        
        // Calculate tides raised on the star by the planet 
        if (stark1 > 0 && starradius > 0 && starmass > 0){
            double startau = fixedstartau; // use particles[0]'s tctl_tau by default 
            // if planet's primary_tau is set, overwrite startau
            double* startauptr = rebx_get_param(rebx, planet->ap, "tctl_primary_tau");
            if (startauptr){
                startau = *startauptr;
            }
            rebx_calculate_tides(planet, star, G, stark1, startau, starOmega);
        }
    }
}

// Calculate potential of conservative piece of tidal interaction
static double rebx_calculate_tides_potential(struct reb_particle* source, struct reb_particle* target, const double G, const double k1){
    const double ms = source->m;
    const double mt = target->m;
    const double Rt = target->r;

    const double mratio = ms/mt; // have already checked for 0 and inf
    const double fac = mratio*k1*Rt*Rt*Rt*Rt*Rt; 
    
    const double dx = target->x - source->x; 
    const double dy = target->y - source->y;
    const double dz = target->z - source->z;
    const double dr2 = dx*dx + dy*dy + dz*dz; 
    
    return -1./2.*G*ms*mt/(dr2*dr2*dr2)*fac;
}

double rebx_tides_constant_time_lag_potential(struct rebx_extras* const rebx){
    if (rebx->sim == NULL){
        rebx_error(rebx, ""); // rebx_error gives meaningful err
        return 0;
    }
    struct reb_simulation* const sim = rebx->sim;
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    const double G = sim->G;
    double H=0.;

    // For conservative tidal potential, we can do stellar tides on all planets separately
    // We only allow dissipative part (primary tau) of tides raised on primary to vary planet by secondary
    // Calculate tides raised on primary
    struct reb_particle* target = &particles[0];// assumes nearly Keplerian motion around a single primary (particles[0])
    if (target->m == 0){                        // No potential with massless primary
        return 0.;
    }
    double* k1 = rebx_get_param(rebx, target->ap, "tctl_k1");
    if (k1 != NULL && target->r != 0){  // tides on star only nonzero if k1 and finite size are set
        for (int i=1; i<N_real; i++){
            struct reb_particle* source = &particles[i]; // planet raising the tides on the star
            if (source->m == 0){
                continue;
            }
            H += rebx_calculate_tides_potential(source, target, G, *k1);
        }
    }

    // Calculate tides raised on the planets
    struct reb_particle* source = &particles[0]; // Source is always the primary (no planet-planet tides)
    for (int i=1; i<N_real; i++){
        struct reb_particle* target = &particles[i]; 
        double* k1 = rebx_get_param(rebx, target->ap, "tctl_k1");
        if (k1 == NULL || target->r == 0 || target->m == 0){
            continue;
        }
        H += rebx_calculate_tides_potential(source, target, G, *k1);
    }

    return H;
}
