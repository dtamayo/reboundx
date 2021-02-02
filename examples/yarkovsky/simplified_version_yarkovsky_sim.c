//
//  simplified_version_yarkovsky_sim.c
//  
//
//  Created by Noah Ferich on 1/14/21.
//

#include "simplified_version_yarkovsky_sim.h"

#include "simple_yarkovsky_sim.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void yarkovsky_effect(struct reb_simulation* const sim);

double tmax = 50000; //in yrs
double radius = 1000; //radius of asteroid (m)
double body_density = 3000; //density of most of the asteroid (kg/m^3)
double lstar = 3.828e26; //luminosity of sun in watts
double c = 299792458;    //speed of light in m/s
double alph = 1;  //alph constant in Veras Yarkovsky equation
double k = 1; //k constant in Veras Yarkovsky equation
double r_conv = 1.495978707e11; //converts AU to m
double t_conv = 31557600; //converts yrs to seconds

int main(int argc, char* argv[]) {
    
    struct reb_simulation* sim = reb_create_simulation(); //creates simulation
    
    sim->G = 4*M_PI*M_PI;  // use units of AU, yr and solar masses
    sim->dt = .05;         //timestep for simulation in yrs
    sim->integrator = REB_INTEGRATOR_WHFAST; //integrator for sim
    
    //following adds star with mass of Sun to sim
    struct reb_particle star = {0};
    star.m = 1.;
    reb_add(sim, star);
    
    //following variables are the orbital elements of only asteroid in sim
    double m = 0;
    double a = .5;
    double e = 0;
    double inc = 0;
    double Omega = 0;
    double omega = 0;
    double f = 0;
    
    //adds asteroid to the sim
    struct reb_particle asteroid = reb_tools_orbit_to_particle(sim->G, star, m, a, e, inc, Omega, omega, f);
    reb_add(sim,asteroid);
    
    sim->additional_forces = yarkovsky_effect; //adds yarkovsky force to sim
    sim->force_is_velocity_dependent = 1; //tells sim force is velocity dependent

    reb_integrate(sim, tmax); //integrates system for tmax years
    
    struct reb_orbit o= reb_tools_particle_to_orbit(sim->G, sim->particles[1], sim->particles[0]); //o gives orbital parameters for asteroid after sim
    
    double final_a = o.a; //final semi-major axis of asteroid after sim
    
    printf("CHANGE IN SEMI-MAJOR AXIS: %1.30f\n", (final_a-a)); //prints difference between the intitial and final semi-major axes of asteroid
    
    /*int i;
    
    double yark_matrix[3][3] = {{1, 2, 3},{1, 5, 3},{6, 2, 3}}; // Same thing as
    
    double i_vector[3][1] = {{1},{2},{3}};
    
    double direction_matrix[3][1];
    
    printf("%f\n", direction_matrix[0][0]);
    
    
    for (i=0; i<3; i++){
         direction_matrix[i][0] = (yark_matrix[i][0]*i_vector[0][0]) + (yark_matrix[i][1]*i_vector[1][0]) + (yark_matrix[i][2]*i_vector[2][0]);
     }

    
    printf("%f\n", direction_matrix[1][0]);*/
    
}

//following calculates yarkovsky effect for asteroid in the sim
void yarkovsky_effect(struct reb_simulation* const sim){
    
    int i; //variable needed for future iteration loops
    
    double v_conv = r_conv/t_conv; //converts AU/yr to m/s
    double a_conv = v_conv/t_conv; //converts AU/yr^2 to m/s^2
    
    struct reb_particle* const particles = sim->particles; //pointer for the particles in the sim
    
    double asteroid_mass = (4/3)*M_PI*(radius*radius*radius)*body_density; //calculates mass of the asteroid from the body density and radius- assume its spherical
    
    double v_vector[3][1] = {{particles[1].vx*v_conv}, {particles[1].vy*v_conv}, {particles[1].vz*v_conv}}; //vector for velocity of asteroid
    
    double r_vector[3][1] = {{particles[1].x*r_conv}, {particles[1].y*r_conv}, {particles[1].z*r_conv}}; //vector for position of asteroid
    
    double distance = sqrt((r_vector[0][0]*r_vector[0][0])+(r_vector[1][0]*r_vector[1][0])+(r_vector[2][0]*r_vector[2][0])); //distance of asteroid from the star
    
    double rdotv = ((r_vector[0][0]*v_vector[0][0])+(r_vector[1][0]*v_vector[1][0])+(r_vector[2][0]*v_vector[2][0]))/(c*distance); //dot product of position and velocity vectors- the term in the denominator is needed when calculating the i vector
    
    double i_vector[3][1];
    
    //loop calculates i vector using equation from Veras, Higuchi, Ida (2019)
    for (i=0; i<3; i++){
     
        i_vector[i][0] = ((1-rdotv)*(r_vector[i][0]/distance))-(v_vector[i][0]/c);
    }
    
    double yarkovsky_magnitude = (radius*radius*lstar)/(4*asteroid_mass*c*distance*distance); //magnitude of the yarkovsky effect for the asteroid
    
    double yark_matrix[3][3] = {{1, 0, 0},{1/4, 1, 0},{0, 0, 1}}; // Same thing as Q in Veras, Higuchi, Ida (2019)
    
    double direction_matrix[3][1];
    
    //loops calcuates a vector which gives the direction of the acceleration created by the Yarkovsky effect
    for (i=0; i<3; i++){
        direction_matrix[i][0] = (yark_matrix[i][0]*i_vector[0][0]) + (yark_matrix[i][1]*i_vector[1][0]) + (yark_matrix[i][2]*i_vector[2][0]);
    }
    
    double yarkovsky_acceleration[3][1];
    
    for (i=0; i<3; i++){
     
        //final result for acceleration created by the Yarkovsky effect- converts it back into units for the sim
        yarkovsky_acceleration[i][0] = (direction_matrix[i][0]*yarkovsky_magnitude)/a_conv;
        
    }
    
    //adds Yarkovsky aceleration to the asteroid's acceleration in the sim
    particles[1].ax += yarkovsky_acceleration[0][0];
    particles[1].ay += yarkovsky_acceleration[1][0];
    particles[1].az += yarkovsky_acceleration[2][0];

}

