//
//  simple_yarkovsky_sim.c
//  
//
//  Created by Noah Ferich on 12/30/20.
//

#include "simple_yarkovsky_sim.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* sim);
void additional_forces(struct reb_simulation* const sim);


double tmax = 1; //in yrs
double radius = 1000; //radius of asteroid (m)
double rotation_period = 15470.9; //how long it takes the body to rotate (sec)
double body_density = 3000; //density of most of the asteroid (kg/m^3)
double C = 680;    //surface heat capacity of asteroid. 832 Karin is an S-type astroid so it's main component is silicon dioxide which has this specific heat value (J/(kg-K))
double c = 299792458;    //speed of light in m/s
double albedo = .017;    //albedo of asteroid- ASSUMED
double alph = 1;  //alph constant in equation
double stef_boltz = 5.670e-8;     //stefan-boltzmann constant (W/(m^2-K^4))
double emissivity = .9; //Found estimate of emissivity of large grains of silicon dioxide (rough approximation)- ASSUMED
double k = 1;
double v_conv = 1.495978707e11/31557600; //converts AU/s to m/s
double r_conv = 1.495978707e11; //converts AU to m
double t_conv = 31557600; //converts yrs to seconds


int main(int argc, char* argv[]) {
    
    struct reb_simulation* sim = reb_create_simulation();

    
    sim->G = 4*M_PI*M_PI;  // use units of AU, yr and solar masses
    sim->dt = .5;
    sim->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.
    
    struct reb_particle sun = {0};
    sun.m = 1.;
    reb_add(sim, sun);
    
    double m = 0.;
    double a = .5;
    double e = 0;
    double inc = 0;
    double Omega = 0;
    double omega = 0;
    double f = 0;
    
    struct reb_particle asteroid = reb_tools_orbit_to_particle(sim->G, sun, m, a, e, inc, Omega, omega, f);
    reb_add(sim,asteroid);
    
    struct reb_particle* const particles = sim->particles;
    
    struct reb_orbit o= reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    
    sim->additional_forces = additional_forces;
    sim->force_is_velocity_dependent = 1;
    
    printf("Intital semi-major axis: %3.15f AU\n", o.a);
    
    reb_integrate(sim, tmax);
    
    printf("Final semi-major axis: %3.15f AU\n", o.a );
    
    printf("%3.15f\n", o.P);

    printf("%2.5f", particles[1].y);
    
}


void additional_forces(struct reb_simulation* const sim){
    
    int i;
    int j;
    
    struct reb_particle* const particles = sim->particles;
    struct reb_orbit o= reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    
    double lstar = 3.828e26; //in watts
    double K = (300*300)/(body_density*C);    //surface thermal conductivity (W/m-K)-ASSUMED REGOLITH_COVERED SURFACE
     double asteroid_mass = (4/3)*M_PI*(radius*radius*radius)*body_density;
    
    double sx = 0.0872;
    double sy = 0;
    double sz = -0.9962;
    double Smag = sqrt((sx*sx)+ (sy*sy) + (sz*sz));
   
    
    double a_conv= v_conv/t_conv; //converts AU/yr^2 to m/s^2
    
    int unit_matrix[3][3] = {{1, 1, 1},{1, 1, 1},{1, 1, 1}};
    
    double distance = sqrt(((particles[1].x*particles[1].x*r_conv*r_conv))+ ((particles[1].y*particles[1].y*r_conv*r_conv)) + ((particles[1].z*particles[1].z*r_conv*r_conv)));
    
    
    double v_vector[3][1] = {{particles[1].vx*v_conv}, {particles[1].vy*v_conv}, {particles[1].vz*v_conv}};
    
    double r_vector[3][1] = {{particles[1].x*r_conv}, {particles[1].y*r_conv}, {particles[1].z*r_conv}};
    
    double hx = (r_vector[1][0]*v_vector[2][0])-(r_vector[2][0]*v_vector[1][0]);
    double hy = (r_vector[2][0]*v_vector[0][0])-(r_vector[0][0]*v_vector[2][0]);
    double hz = (r_vector[0][0]*v_vector[1][0])-(r_vector[1][0]*v_vector[0][0]);
    double Hmag = sqrt((hx*hx)+ (hy*hy) + (hz*hz));
    
    double R1s[3][3] = {{0, -sz*(1/Smag), sy*(1/Smag)},{sz*(1/Smag), 0, -sx*(1/Smag)},{-sy*(1/Smag), sx*(1/Smag), 0}};
    double R2s[3][3] = {{sx*sx*(1/(Smag*Smag)), sx*sy*(1/(Smag*Smag)), sx*sz*(1/(Smag*Smag))},{sx*sy*(1/(Smag*Smag)), sy*sy*(1/(Smag*Smag)), sy*sz*(1/(Smag*Smag))},{sx*sz*(1/(Smag*Smag)), sy*sz*(1/(Smag*Smag)), sz*sz*(1/(Smag*Smag))}};
    double R1h[3][3] = {{0, -hz*(1/Hmag), hy*(1/Hmag)},{hz*(1/Hmag), 0, -hx*(1/Hmag)},{-hy*(1/Hmag), hx*(1/Hmag), 0}};
    double R2h[3][3] = {{hx*hx*(1/(Hmag*Hmag)), hx*hy*(1/(Hmag*Hmag)), hx*hz*(1/(Hmag*Hmag))},{hx*hy*(1/(Hmag*Hmag)), hy*hy*(1/(Hmag*Hmag)), hy*hz*(1/(Hmag*Hmag))},{hx*hz*(1/(Hmag*Hmag)), hy*hz*(1/(Hmag*Hmag)), hz*hz*(1/(Hmag*Hmag))}};
    
    double rdotv = ((r_vector[0][0]*v_vector[0][0])+(r_vector[1][0]*v_vector[1][0])+(r_vector[2][0]*v_vector[2][0]))/(c*distance);
    
    double i_vector[3][1];
    
    for (i=0; i<3; i++){
     
        i_vector[i][0] = (((1-rdotv)/distance)*r_vector[i][0])-(v_vector[i][0]/c);
    }
    
    double tanPhi = 1/(1+(.5*pow((stef_boltz*emissivity)/(M_PI*M_PI*M_PI*M_PI*M_PI), .25))*sqrt(rotation_period/(C*K*body_density))*pow((lstar*alph)/(distance*distance), .75));
    
    double tanEpsilon = 1/(1+(.5*pow((stef_boltz*emissivity)/(M_PI*M_PI*M_PI*M_PI*M_PI), .25))*sqrt((o.P*t_conv)/(C*K*body_density))*pow((lstar*alph)/(distance*distance), .75));
    
    double Phi = atan(tanPhi);
    double Epsilon = atan(tanEpsilon);
    
    double Rys[3][3];
    
    for (i=0, j=0; i<3 && j<3; i++, j++){
        Rys[i][j] = (cos(Phi)*unit_matrix[i][j]) + sin(Phi)*R1s[i][j] + (1-cos(Phi))*R2s[i][j];
    }
    
    double Ryh[3][3];
    
    for (i=0, j=0; i<3 && j<3; i++, j++){
        Ryh[i][j] = (cos(Epsilon)*unit_matrix[i][j]) + sin(Epsilon)*R1h[i][j] + (1-cos(Epsilon))*R2h[i][j];
    }
    
    double yarkovsky_magnitude = ((radius*radius)*lstar*alph)/(4*asteroid_mass*c*(distance*distance));
    
    //printf("%1.20f\n", yarkovsky_magnitude);
    
    double yark_matrix[3][3];
    
    for (i=0, j=0; i<3 && j<3; i++, j++){
            yark_matrix[i][j] = Rys[i][0]*Ryh[0][j] + Rys[i][1]*Ryh[1][j] + Rys[i][2]*Ryh[2][j];
        }
    
    double direction_matrix[3][1];
    
    for (i=0; i<3; i++){
        direction_matrix[i][0] = yark_matrix[i][0]*i_vector[0][0] + yark_matrix[i][1]*i_vector[1][0] + yark_matrix[i][2]*i_vector[2][0];
    }
    
    printf("%1.7f\n",direction_matrix[0][0]);
    printf("%1.7f\n",direction_matrix[1][0]);
    printf("%1.7f\n",direction_matrix[2][0]);
    
    double yarkovsky_acceleration[3][1];
    
    for (i=0; i<3; i++){
     
        yarkovsky_acceleration[i][0] = (direction_matrix[i][0]*yarkovsky_magnitude)/a_conv;
        
    }

    
    particles[1].ax += yarkovsky_acceleration[0][0];
    particles[1].ay += yarkovsky_acceleration[1][0];
    particles[1].az += yarkovsky_acceleration[2][0];
    

}

