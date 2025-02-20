/**
 * @file    tides_dynamical.c
 * @brief   Update body's orbital and modal evolution due to the presence of dynamical tides.
 * @author  Donald J. Liveoak <donaldliveoak1@gmail.com>
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
 * Authors                 D. Liveoak & S. Millholland
 * Implementation Paper    
 * Based on                `Vick et al. 2019 <https://academic.oup.com/mnras/article/484/4/5645/5306464/>`_. 
 * C Example               
 * Python Example          
 * ======================= ===============================================
 * 
 * This updates body's orbital and modal evolution due to the presence of dynamical tides.
 * Particles are modeled by a gamma=4/3 polytrope, and the f-mode is evolved at each pericentre passage.
 * The dissipation of orbital energy due to dynamical tides is modeled as an angular momentum-conserving kick at periapse.
 * When mode energy grows to exceed `td_E_max`, it is non-linearly dissipated in one orbital period to `td_E_resid`.
 * To isolate the effects of chaotic model evolution, one can set `dP_hat_crit` to disable dynamical tides whenever chaos is unlikely (see Vick et al. (2019))
 * 
 * 
 * **Particle Parameters**
 *
 * One can pick and choose which particles have which parameters set.  
 * For each particle, any unset parameter is replaced by its default value.
 * Particles with index 0 will not experience dynamical tides.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

struct rebx_tides_dynamical_params rebx_calculate_tides_dynamical_params(struct reb_simulation* const sim, struct reb_particle* p, struct reb_particle* primary)
{

    struct rebx_extras* const rebx = sim->extras;
    const double EulerConstant = 2.718281828459;

    // Calculate orbital elements
    struct reb_orbit o = reb_orbit_from_particle(sim->G, *p, *primary);
    double e = o.e;
    double a = o.a;
    double n = o.n;
    double P = o.P;

    // Calculate some useful distances
    double R = p->r; // radius of planet
    double R_tide = R * pow(primary->m / p->m, 0.33333); // tidal radius
    double R_p = a * (1 - e); // pericenter distance
    double eta = R_p / R_tide; // pericenter distance in units of tidal radius


    if (eta <= 2.5) // Tidal disruption occured
    {
        struct rebx_tides_dynamical_params toReturn;
        toReturn.dP = 0;
        toReturn.dE_alpha = 0;
        toReturn.sigma = 0;

        return toReturn;
    }

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

    struct rebx_tides_dynamical_params toReturn;
    toReturn.dP = dP;
    toReturn.dE_alpha = dE_alpha;
    toReturn.sigma = sigma;

    return toReturn;
}

struct rebx_tides_dynamical_mode rebx_calculate_tides_dynamical_mode_evolution(double old_real, double old_imag, double dc_tilde, double P, double sigma)
{
    double new_real = (old_real + dc_tilde) * cos(sigma * P) + old_imag * sin(sigma * P);
    double new_imag = -(old_real + dc_tilde) * sin(sigma * P) + old_imag * cos(sigma * P);

    struct rebx_tides_dynamical_mode mode;
    mode.real = new_real;
    mode.imag = new_imag;

    return mode;
}

double rebx_calculate_tides_dynamical_drag_integral(double e, double n)
{
    if (n == 4)
    {
        return M_PI * (2 + 7 * e*e + e*e*e*e) / 2;
    }
    if (n == 10)
    {
        return M_PI * (128 + 2944*e*e + 10528*pow(e, 4) + 8960*pow(e, 6) + 1715*pow(e, 8) + 35*pow(e, 10)) / 128;
    }
    return 0;
}

void rebx_tides_dynamical(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N)
{
   	struct rebx_extras* const rebx = sim->extras;
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
    if (rebx_get_param(rebx, p->ap, "td_dP_hat") == NULL)
    {
        rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_dP_hat", 0);
    }
    if (rebx_get_param(rebx, p->ap, "td_drag_coef") == NULL)
    {
        rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_drag_coef", 0);
    }
    if (rebx_get_param(rebx, p->ap, "td_drag_exp") == NULL)
    {
        rebx_set_param_int(rebx, (struct rebx_node**)&p->ap, "td_drag_exp", 4);
    }
    if (rebx_get_param(rebx, p->ap, "td_last_periapse") == NULL)
    {
        rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_last_periapse", 0);
    }
    if (rebx_get_param(rebx, p->ap, "td_debug_Eb_last") == NULL)
    {
        double EB0 = -sim->G * p->m * source->m / (2 * o.a);
        rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_debug_Eb_last", EB0);
    }
    if (rebx_get_param(rebx, p->ap, "td_dc_tilde_last") == NULL)
    {
        rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_dc_tilde_last", 0);
    }

    if (rebx_get_param(rebx, p->ap, "td_M_last") != NULL)
    {       
        double* M_last = rebx_get_param(rebx, p->ap, "td_M_last");
        double drag = 0;
        double* last_periapse_time = rebx_get_param(rebx, p->ap, "td_last_periapse");
        if ((o.M >= M_PI && *M_last < M_PI) && o.M - M_PI <= 1 && sim->t - *last_periapse_time >= sim->dt)        //if (o.M < *M_last)
        {
            // Count periapse passages
            int* num_periapse = rebx_get_param(rebx, p->ap, "td_num_periapse");
            rebx_set_param_int(rebx, (struct rebx_node**)&p->ap, "td_num_periapse", *num_periapse + 1);

            double* dP_crit = rebx_get_param(rebx, p->ap, "td_dP_crit");
            struct rebx_tides_dynamical_params dynamical_params = rebx_calculate_tides_dynamical_params(sim, p, source);
            double dP = dynamical_params.dP;
            double dE_alpha = dynamical_params.dE_alpha;
            
            rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_dP_hat", dP);
            rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_dE_last", dE_alpha);

            // If system is in chaotic regime, evolve dynamical tides
            if (dP >= *dP_crit)
            {
                // Calculate map parameters
                double* EB0 = rebx_get_param(rebx, p->ap, "td_EB0");
                double EBk = -sim->G * p->m * source->m / (2 * o.a);
                double dc_tilde = pow(dE_alpha / -*EB0, 0.5);
                double dE_alpha_tilde = dE_alpha / -*EB0;
                double* c_real = rebx_get_param(rebx, p->ap, "td_c_real"); 
                double* c_imag = rebx_get_param(rebx, p->ap, "td_c_imag");
                rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_dc_tilde_last", dc_tilde);

                // Evolve modes
                double sigma = dynamical_params.sigma;
                struct rebx_tides_dynamical_mode new_modes = rebx_calculate_tides_dynamical_mode_evolution(*c_real, *c_imag, dc_tilde, o.P, sigma);
                double dEb = (-*EB0) * (new_modes.real*new_modes.real + new_modes.imag*new_modes.imag - *c_real* *c_real - *c_imag * *c_imag);
                double* EB_last = rebx_get_param(rebx, p->ap, "td_debug_Eb_last");
                rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_debug_Eb_last", EBk);

                double* E_max = rebx_get_param(rebx, p->ap, "td_E_max"); 
                double* E_resid = rebx_get_param(rebx, p->ap, "td_E_resid");
                if (-(pow(new_modes.real, 2) + pow(new_modes.imag, 2)) * *EB0 >= *E_max)
                {
                    double E_dis_ratio = -*E_resid / *EB0;
                    new_modes.real = pow(E_dis_ratio / (1 + pow(new_modes.imag, 2) / pow(new_modes.real, 2)), 0.5);
                    new_modes.imag = pow(E_dis_ratio / (1 + pow(new_modes.real, 2) / pow(new_modes.imag, 2)), 0.5);

                }

                rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_c_real", new_modes.real);
                rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_c_imag", new_modes.imag);  
                rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_last_periapse", sim->t);


                // Compute drag parameter
                int* n = rebx_get_param(rebx, p->ap, "td_drag_exp");
                double I = rebx_calculate_tides_dynamical_drag_integral(o.e, *n);
                
                drag = dEb * pow((o.a) * (1 - o.e * o.e), *n - 0.5) / (2 * pow(sim->G * (p->m + source->m), 0.5) * I);
            }
            rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_drag_coef", drag);
        }
    }

    rebx_set_param_double(rebx, (struct rebx_node**)&p->ap, "td_M_last", o.M);    

    double* drag_coef = rebx_get_param(rebx, p->ap, "td_drag_coef");
    int* drag_exp = rebx_get_param(rebx, p->ap, "td_drag_exp");

    // Compute CoM
    double comx = 0;
    double comy = 0;
    double comz = 0;
    double comvx = 0;
    double comvy = 0;
    double comvz = 0;
    double total_m = 0;

    for (int i = 0; i < 2; i++)
    {
        comx += particles[i].m * particles[i].x;
        comy += particles[i].m * particles[i].y;
        comz += particles[i].m * particles[i].z;
        comvx += particles[i].m * particles[i].vx;
        comvy += particles[i].m * particles[i].vy;
        comvz += particles[i].m * particles[i].vz;
        total_m += particles[i].m;
    }

    double x = particles[1].x - comx / total_m;
    double y = particles[1].y - comy / total_m;
    double z = particles[1].z - comz / total_m;
    double vx = particles[1].vx - comvx / total_m;
    double vy = particles[1].vy - comvy / total_m;
    double vz = particles[1].vz - comvz / total_m;

    double r = pow(x*x + y*y + z*z, 0.5);
    double Fx = -*drag_coef * vx / pow(r, *drag_exp);
    double Fy = -*drag_coef * vy / pow(r, *drag_exp);
    double Fz = -*drag_coef * vz / pow(r, *drag_exp);

    // Apply drag to particle
    particles[1].ax += Fx / particles[1].m;
    particles[1].ay += Fy / particles[1].m;
    particles[1].az += Fz / particles[1].m;

    // Apply equal and opposite force to primary
    particles[0].ax -= Fx / particles[0].m;
    particles[0].ay -= Fy / particles[0].m;
    particles[0].az -= Fz / particles[0].m;
}
