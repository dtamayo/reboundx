/**
 * @file    modify_orbits_direct.c
 * @brief   Update orbital with prescribed timescales by directly changing orbital elements after each timestep.
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
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
 * $Orbit Modifications$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_. 
 * Based on                `Lee & Peale 2002 <http://labs.adsabs.harvard.edu/adsabs/abs/2002ApJ...567..596L/>`_. 
 * C Example               :ref:`c_example_modify_orbits`
 * Python Example          `Migration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Migration.ipynb>`_,
 *                         `EccAndIncDamping.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/EccAndIncDamping.ipynb>`_.
 * ======================= ===============================================
 * 
 * This updates particles' positions and velocities between timesteps to achieve the desired changes to the osculating orbital elements (exponential growth/decay for a, e, inc, linear progression/regression for Omega/omega.
 * This nicely isolates changes to particular osculating elements, making it easier to interpret the resulting dynamics.  
 * One can also adjust the coupling parameter `p` between eccentricity and semimajor axis evolution, as well as whether the damping is done on Jacobi, barycentric or heliocentric elements.
 * Since this method changes osculating (i.e., two-body) elements, it can give unphysical results in highly perturbed systems.
 * 
 * **Effect Parameters**
 *
 * If p is not set, it defaults to 0.  If coordinates not set, defaults to using Jacobi coordinates.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * p (double)                   No          Coupling parameter between eccentricity and semimajor axis evolution
 *                                          (see Deck & Batygin 2015). `p=0` corresponds to no coupling, `p=1` to
 *                                          eccentricity evolution at constant angular momentum.
 * coordinates (enum)           No          Type of elements to use for modification (Jacobi, barycentric or particle).
 *                                          See the examples for usage.
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * One can pick and choose which particles have which parameters set.  
 * For each particle, any unset parameter is ignored.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * tau_a (double)               No          Semimajor axis exponential growth/damping timescale
 * tau_e (double)               No          Eccentricity exponential growth/damping timescale
 * tau_inc (double)             No          Inclination axis exponential growth/damping timescale
 * tau_Omega (double)           No          Period of linear nodal precession/regression
 * tau_omega (double)           No          Period of linear apsidal precession/regression
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

double[] rebx_calculate_tides_dynamical_params(struct reb_simulation* const sim, struct reb_particle* p, struct reb_particle* primary)
{

    struct rebx_extras* const rebx = sim->extras;
    const double EulerConstant = 2.718281828459;

    // Calculate orbit
    struct reb_orbit o = reb_orbit_from_particle(sim->G, *p, *primary);
    double e = o.e;
    double a = o.a;
    double n = o.n;
    double P = o.P;

    // Calculate some useful distances
    double R = p->r; // radius of planet
    double R_tide = R * pow(primary->m / p->m, 1 / 3); // tidal radius
    double R_p = a * (1 - e); // pericenter distance
    double eta = R_p / R_tide; // pericenter distance in units of tidal radius

    // Timescales/frequencies
    double Omega_peri = pow(sim->G * (p->m + primary->m) / (R_p * R_p * R_p), 0.5); // pericenter frequency
    double time_unit = pow(sim->G * p->m / (R * R * R), 0.5); // default units for mode parameters

    // Calculate pseudo-synchronous orbital frequency
    double f2 = 1 + (15 / 2) * pow(e, 2) + (45 / 8) * pow(e, 4) + (5 / 16) * pow(e, 6);
    double f5 = 1 + 3 * pow(e, 2) + (3 / 8) * pow(e, 4);
    double Omega_s = n * f2 / (pow(1 - e * e, 1.5) * f5);

    // Calculate f-mode parameters, gamma=2 polytrope (see Vick et al. (2019))
    double omega = (1.22 - Omega_s / time_unit) * time_unit;
    double sigma = (1.22 + Omega_s / time_unit) * time_unit;
    double epsilon = 1.22 * time_unit;
    double Q = 0.56; // overlap integral

    // Calculate K_22 and T
    double z = pow(2, 0.5) * sigma / Omega_peri;
    double K_22 = 2 * pow(z, 1.5) * pow(eta, 1.5) * pow(EulerConstant, -2 * z / 3) * (1 - pow(M_PI, 0.5) / (4 * pow(z, 0.5))) / (pow(15, 0.5));
    double T = 2 * M_PI * M_PI * Q * Q * K_22 * K_22 * sigma / epsilon;

    // Calculate change in mode energy, assuming 0 mode amplitude
    double dE_alpha = sim->G * primary->m * primary->m * pow(R, 5) * T / pow(R_p, 6);

    // Calculate dP
    double* EB0 = rebx_get_param(rebx, p->ap, "td_EB0");
    double* c_real = rebx_get_param(rebx, p->ap, "td_c_real"); 
    double* c_imag = rebx_get_param(rebx, p->ap, "td_c_imag");
    double maxE = dE_alpha + 2 * pow(-dE_alpha * (pow(*c_real, 2) + pow(*c_imag, 2)) * *EB0, 0.5);
    double EBk = -sim->G * p->m * primary->m / (2 * a);
    double dP = 1.5 * sigma * P * maxE / (-EBk);

    return {dP, dE_alpha, sigma};
}

double rebx_calculate_tides_dynamical_calc_mode_evolution(double[] modes, double dc_tilde, double P, double sigma)
{
    double old_real = modes[0];
    double old_imag = modes[1];
    double new_real = (modes[0] + dc_tilde) * cos(sigma * P) + modes[1] * sin(sigma * P);
    double new_imag = -(modes[0] + dc_tilde) * sin(sigma * P) + modes[1] * cos(sigma * P);
    return {new_real, new_imag};
}


void rebx_tides_dynamical(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){

   	struct rebx_extras* const rebx = sim->extras;

    // update the momentum of the first particle
    struct reb_particle* const source = &sim->particles[0];
    struct reb_particle* const p = &sim->particles[1];

    struct reb_orbit o = reb_orbit_from_particle(sim->G, *p, *source);

    // Set default parameter values
    if (rebx_get_param(rebx, p->ap, "td_EB0") == NULL)
    {
        double EB0 = -sim->G * p->m * source->m / (2 * o.a);
        rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_EB0", EB0);    
    }
    if (rebx_get_param(rebx, p->ap, "td_num_periapse") == NULL)
    {
        rebx_set_param_int(rebx, (struct rebx_node**)&p->ap, "td_num_periapse", 0);    
    }
    if (rebx_get_param(rebx, p->ap, "td_c_real") == NULL)
    {
        rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_c_real", 0);    
    }
    if (rebx_get_param(rebx, p->ap, "td_c_imag") == NULL)
    {
        rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_c_imag", 0);    
    }
    if (rebx_get_param(rebx, p->ap, "td_dP_crit") == NULL)
    {
        rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_dP_crit", 0.01);    
    }
    if (rebx_get_param(rebx, p->ap, "td_E_max") == NULL)
    {
        double E_bind = sim->G * p->m * p->m / p->r;
        rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_E_max", E_bind / 10);    
    }
    if (rebx_get_param(rebx, p->ap, "td_E_resid") == NULL)
    {
        double E_bind = sim->G * p->m * p->m / p->r;
        rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_E_resid", E_bind / 1000);    
    }


    if (rebx_get_param(rebx, p->ap, "td_M_last") != NULL)
    {       
        double* M_last = rebx_get_param(rebx, p->ap, "td_M_last");

        // Periapse detection
        if (o.M < *M_last)
        {
            // Count periapse passages
            int* num_periapse = rebx_get_param(rebx, p->ap, "td_num_periapse");
            rebx_set_param_int(rebx, (struct rebx_node**)&p->ap, "td_num_periapse", *num_periapse + 1);

            double* dP_crit = rebx_get_param(rebx, p->ap, "dP_crit");
            double[] dynamical_params = rebx_calculate_tides_dynamical_params(sim, p, source);
            double dP = dynamical_parmas[0];
            rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_dP_hat", dP);

            // If system is in chaotic regime, evolve dynamical tides
            if (dP >= dP_crit)
            {
                double* EB0 = rebx_get_param(rebx, p->ap, "td_EB0");
                double EBk = -sim->G * p->m * source->m / (2 * o.a);
                double dE_alpha = dynamical_params[1];
                double dc_tilde = pow(dE_alpha / -EB0, 0.5);
                double dE_alpha_tilde = dE_alpha / -EB0;
                double* c_real = rebx_get_param(rebx, p->ap, "td_c_real"); 
                double* c_imag = rebx_get_param(rebx, p->ap, "td_c_imag");

                // Calculate new orbital energy
                double EB_new = EBk - (-Ebk) * (dE_alpha_tilde + 2 * pow(dE_alpha_tilde, 0.5) * c_real);
                double E_ratio = EBk / E_B;

                // If amplitude is sufficiently high, non-linear dissipation
                double* E_max = rebx_get_param(rebx, p->ap, "td_E_max"); 
                double* E_resid = rebx_get_param(rebx, p->ap, "td_E_resid");
                if (-(pow(*c_real, 2) + pow(*c_imag, 2)) * *EB0 >= *E_max)
                {
                    double E_dis_ratio = -E_resid / *EB0;
                    *c_real = pow(E_dis_ratio / (1 + pow(*c_imag, 2) / pow(*c_real, 2)));
                    *c_imag = pow(E_dis_ratio / (1 + pow(*c_real, 2) / pow(*c_imag, 2)));
                }

                // Calculate new orbital elements
                double a_prime = o.a * E_ratio;
                double e_prime = pow(1 - (1 / E_ratio) * (1 - o.e * o.e), 0.5);
                double P_prime = o.P * pow(E_ratio, -1.5);

                // Evolve modes
                double[] new_modes = rebx_calculate_tides_dynamical_calc_mode_evolution({*c_real, *c_imag}, dc_tilde, P, sigma);
                rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_c_real", new_modes[0]);  
                rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_c_imag", new_modes[1]);  
            }
        }

    }

    
    rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_M_last", o.M);    
    
}
