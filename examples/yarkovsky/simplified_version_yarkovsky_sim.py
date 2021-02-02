import rebound
import numpy as np

#Simulation begins here
sim = rebound.Simulation()

sp = sim.particles #simplifies way to access particles parameters 

sim.units = ('yr', 'AU', 'Msun') #changes simulation and G to units of solar masses, years, and AU  
sim.integrator = "whfast" #integrator for sim
sim.dt = .05 #timestep for sim


#Adds Sun and outer planets to simulation
sim.add(m=1) 
sim.add(a=.5, f=0, Omega=0, omega=0, e=0, inc=0, m=0)

#Moves all particles to center of momentum frame
sim.move_to_com()

#Gives orbital information before the simulation begins
print("\n***INITIAL ORBITS:***")
for orbit in sim.calculate_orbits():
    print(orbit)

tmax=50000 # in yrs
radius = 1000   #radius of asteroid (m)
lstar = 3.828e26 #luminosity of sun in watts
body_density = 3000   #density of asteroid (kg/m^3)
c = 299792458   #speed of light (m/s)  
alph = 1  #alph constant in equation
mass = ((4/3)*np.pi*(radius**3))*body_density #calculates mass of asteroid in kg
k = 1 #k constant in Veras Yarkovsky equation
r_conv = 1.495978707e11 #coverts AU to m
t_conv = 31557600 #converts yrs to sec
v_conv = r_conv/t_conv #converts AU/yr to m/s
a_conv = v_conv/t_conv #converts AU/yr^2 to m/s^2

#following calculates yarkovsky effect for asteroid in the sim
def yarkovsky_effect(reb_sim):
    
    
    v_vector = np.array([[(sp[1].vx*v_conv)], [(sp[1].vy*v_conv)], [(sp[1].vz*v_conv)]]) #vector for velocity of asteroid

    r_vector = np.array([[(sp[1].x*r_conv)], [(sp[1].y*r_conv)], [(sp[1].z*r_conv)]]) #vector for position of asteroid

    
    distance = ((r_vector[0]**2)+(r_vector[1]**2)+(r_vector[2])**2)**(1/2) #distance of asteroid from the star


    rdotv = ((r_vector[0][0]*v_vector[0][0])+(r_vector[1][0]*v_vector[1][0])+(r_vector[2][0]*v_vector[2][0]))/(c*distance) #dot product of position and velocity vectors- the term in the denominator is needed when calculating the i vector

    i_vector = ((1-rdotv)*(r_vector/distance))-(v_vector/c) #calculates i vector using equation from Veras, Higuchi, Ida (2019)

    yarkovsky_magnitude = ((radius**2)*lstar)/(4*mass*c*(distance**2)) #magnitude of the yarkovsky effect for the asteroid

    Yark_matrix =  np.array([[1, 0, 0], [1/4, 1, 0 ], [0, 0, 1]]) #Same thing as Q in Veras, Higuchi, Ida (2019)
    
    Direction_matrix = Yark_matrix.dot(i_vector) #vector which gives the direction of the acceleration created by the Yarkovsky effect
    
    yarkovsky_acceleration = (yarkovsky_magnitude*Direction_matrix)/a_conv #final result for acceleration created by the Yarkovsky effect- converts it back into units for the sim
    

    #adds Yarkovsky aceleration to the asteroid's acceleration in the sim
    sp[1].ax += yarkovsky_acceleration[0][0]				
    sp[1].ay += yarkovsky_acceleration[1][0]
    sp[1].az += yarkovsky_acceleration[2][0]

sim.additional_forces = yarkovsky_effect  # adds force to simulation 
sim.force_is_velocity_dependent = 1 # tells sim the force depends on particle's velocity

a = .5 #starting semi-major axis for the asteroid
            
sim.integrate(tmax) #integrates system for tmax years       
 
final_a = sp[1].a #semi-major axis of asteroid after the sim    
                      
print("CHANGE IN SEMI-MAJOR AXIS", final_a-a) #prints difference between the intitial and final semi-major axes of asteroid
