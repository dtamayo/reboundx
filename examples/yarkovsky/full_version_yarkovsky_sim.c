//
//  simple_yarkovsky_sim.c
//  
//
//  Created by Noah Ferich on 12/30/20.
//


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* sim);
void additional_forces(struct reb_simulation* const sim);


double tmax = 10; //in yrs
double radius = 1000; //radius of asteroid (m)
    double lstar = 3.828e26; //in watts
double rotation_period = 15470.9; //how long it takes the body to rotate (sec)
double body_density = 3000; //density of most of the asteroid (kg/m^3)
double C = 680;    //surface heat capacity of asteroid. 832 Karin is an S-type astroid so it's main component is silicon dioxide which has this specific heat value (J/(kg-K))
double c = 299792458;    //speed of light in m/s
double albedo = .017;    //albedo of asteroid- ASSUMED
double stef_boltz = 5.670e-8;     //stefan-boltzmann constant (W/(m^2-K^4))
double emissivity = .9; //Found estimate of emissivity of large grains of silicon dioxide (rough approximation)- ASSUMED
double k = 1;
double r_conv = 1.495978707e11; //converts AU to m
double t_conv = 31557600; //converts yrs to seconds


int main(int argc, char* argv[]) {
    
    struct reb_simulation* sim = reb_create_simulation();

    
    sim->G = 4*M_PI*M_PI;  // use units of AU, yr and solar masses
    sim->dt = .1;
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
    
    sim->additional_forces = additional_forces;
    sim->force_is_velocity_dependent = 1;
    
    reb_integrate(sim, tmax);
    
    struct reb_orbit o= reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]); //o gives orbital parameters for asteroid after sim
    
    double final_a = o.a; //final semi-major axis of asteroid after sim
    
    printf("CHANGE IN SEMI-MAJOR AXIS: %1.30f\n", (final_a-a)); //prints difference between the intitial and final semi-major axes of asteroid
    
}


void additional_forces(struct reb_simulation* const sim){
    
    int i;
    int j;
    
    double v_conv = r_conv/t_conv; //converts AU/yr to m/s
    double a_conv = v_conv/t_conv; //converts AU/yr^2 to m/s^2
    
    struct reb_particle* const particles = sim->particles;
    struct reb_orbit o= reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]);
    
    double K = (300.0*300.0)/(body_density*C);    //surface thermal conductivity (W/m-K)-ASSUMED REGOLITH_COVERED SURFACE
    double asteroid_mass = (4*M_PI*(radius*radius*radius)*body_density)/3;
    double q_yar = 1-albedo;
    
    double sx = 0.0872;
    double sy = 0.0;
    double sz = -0.9962;
    double Smag = sqrt((sx*sx)+ (sy*sy) + (sz*sz));
   
    
    double unit_matrix[3][3] = {{1.0, 1.0, 1.0},{1.0, 1.0, 1.0},{1.0, 1.0, 1.0}};
    
    double v_vector[3][1] = {{particles[1].vx*v_conv}, {particles[1].vy*v_conv}, {particles[1].vz*v_conv}}; //vector for velocity of asteroid
    
    double r_vector[3][1] = {{particles[1].x*r_conv}, {particles[1].y*r_conv}, {particles[1].z*r_conv}}; //vector for position of asteroid
    
    double distance = sqrt((r_vector[0][0]*r_vector[0][0])+(r_vector[1][0]*r_vector[1][0])+(r_vector[2][0]*r_vector[2][0])); //distance of asteroid from the star
    
    
    double hx = (r_vector[1][0]*v_vector[2][0])-(r_vector[2][0]*v_vector[1][0]);
    double hy = (r_vector[2][0]*v_vector[0][0])-(r_vector[0][0]*v_vector[2][0]);
    double hz = (r_vector[0][0]*v_vector[1][0])-(r_vector[1][0]*v_vector[0][0]);
    double Hmag = sqrt((hx*hx)+ (hy*hy) + (hz*hz));
    
    double inv_smag = 1.0/Smag;
    double inv_mag_sqrd = 1.0/(Smag*Smag);
    double inv_hmag = 1.0/Hmag;
    double inv_hmag_sqrd = 1.0/(Hmag*Hmag);
    
    double R1s[3][3] = {{0.0, -sz*inv_smag, sy*inv_smag},{sz*inv_smag, 0.0, -sx*inv_smag},{-sy*inv_smag, sx*inv_smag, 0.0}};
    double R2s[3][3] = {{sx*sx*inv_mag_sqrd, sx*sy*inv_mag_sqrd, sx*sz*inv_mag_sqrd},{sx*sy*inv_mag_sqrd, sy*sy*inv_mag_sqrd, sy*sz*inv_mag_sqrd},{sx*sz*inv_mag_sqrd, sy*sz*inv_mag_sqrd, sz*sz*inv_mag_sqrd}};
    double R1h[3][3] = {{0.0, -hz*inv_hmag, hy*inv_hmag},{hz*inv_hmag, 0.0, -hx*inv_hmag},{-hy*inv_hmag, hx*inv_hmag, 0.0}};
    double R2h[3][3] = {{hx*hx*inv_hmag_sqrd, hx*hy*inv_hmag_sqrd, hx*hz*inv_hmag_sqrd},{hx*hy*inv_hmag_sqrd, hy*hy*inv_hmag_sqrd, hy*hz*inv_hmag_sqrd},{hx*hz*inv_hmag_sqrd, hy*hz*inv_hmag_sqrd, hz*hz*inv_hmag_sqrd}};
    
    
    double rdotv = ((r_vector[0][0]*v_vector[0][0])+(r_vector[1][0]*v_vector[1][0])+(r_vector[2][0]*v_vector[2][0]))/(c*distance); //dot product of position and velocity vectors- the term in the denominator is needed when calculating the i vector
    
    double i_vector[3][1];
    
    //loop calculates i vector using equation from Veras, Higuchi, Ida (2019)
    for (i=0; i<3; i++){
     
        i_vector[i][0] = ((1-rdotv)*(r_vector[i][0]/distance))-(v_vector[i][0]/c);
    }
    
    double tanPhi = 1.0/(1.0+(.5*pow((stef_boltz*emissivity)/(M_PI*M_PI*M_PI*M_PI*M_PI), .25))*sqrt(rotation_period/(C*K*body_density))*pow((lstar*q_yar)/(distance*distance), .75));
    
    double tanEpsilon = 1.0/(1.0+(.5*pow((stef_boltz*emissivity)/(M_PI*M_PI*M_PI*M_PI*M_PI), .25))*sqrt((o.P*t_conv)/(C*K*body_density))*pow((lstar*q_yar)/(distance*distance), .75));
    
    double Phi = atan(tanPhi);
    double Epsilon = atan(tanEpsilon);
    
    double cos_phi = cos(Phi);
    double sin_phi = sin(Phi);
    double sin_epsilon = sin(Epsilon);
    double cos_epsilon = cos(Epsilon);
    
    double Rys[3][3];
    
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            Rys[i][j] = (cos_phi*unit_matrix[i][j]) + (sin_phi*R1s[i][j]) + ((1.0-cos_phi)*R2s[i][j]);
        }
    }
    
    
    
    double Ryh[3][3];
    
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            Ryh[i][j] = (cos_epsilon*unit_matrix[i][j]) - (sin_epsilon*R1h[i][j]) + ((1-cos_epsilon)*R2h[i][j]);
        }
    }
    
    
    double yarkovsky_magnitude = (radius*radius*lstar*q_yar)/(4*asteroid_mass*c*distance*distance); //magnitude of the yarkovsky effect for the asteroid
    
    
    double yark_matrix[3][3];
    
    for (i=0; i<3; i++){
        for(j=0; j<3; j++){
            yark_matrix[i][j] = (Ryh[i][0]*Rys[0][j]) + (Ryh[i][1]*Rys[1][j]) + (Ryh[i][2]*Rys[2][j]);
            }
        }
    
    double direction_matrix[3][1];
    
    for (i=0; i<3; i++){
    direction_matrix[i][0] = (yark_matrix[i][0]*i_vector[0][0]) + (yark_matrix[i][1]*i_vector[1][0]) + (yark_matrix[i][2]*i_vector[2][0]);
    }
    
    
    double yarkovsky_acceleration[3][1];
    
    for (i=0; i<3; i++){
     
        yarkovsky_acceleration[i][0] = (direction_matrix[i][0]*yarkovsky_magnitude)/a_conv;
        
    }
    
    particles[1].ax += yarkovsky_acceleration[0][0];
    particles[1].ay += yarkovsky_acceleration[1][0];
    particles[1].az += yarkovsky_acceleration[2][0];
    

}

