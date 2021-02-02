import mesa_reader as mr
import time
import rebound
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import os


sim = rebound.Simulation()
sp = sim.particles

#Changes simulation and G to units of solar masses, years, and AU
sim.units = ('yr', 'AU', 'Msun') 
sim.add(m=1, x=0, y=0, z=0, hash = 'sun')
sim.add(a=.5)
sim.dt = .01



#Most if not all of the following will be in kms units

radius = 1000     #radius of asteroid (m)
rotation_period = 15470.9    #how long it takes body to rotate (sec)
C = 680    #surface heat capacity of asteroid. 832 Karin is an S-type astroid so it's main component is silicon dioxide which has this specific heat value (J/(kg-K))
lsun = 3.828e26
body_density = 3000   #density of most of asteroid (kg/m^3)
K = (300**2)/(body_density*C)    #surface thermal conductivity (W/m-K)-ASSUMED REGOLITH_COVERED SURFACE
c = 299792458    #speed of light   
albedo = .017    #albedo of asteroid- ASSUMED
alph = 1  #alph constant in equation
stef_boltz = 5.670e-8     #stefan-boltzmann constant (W/(m^2-K^4)) 
emissivity = .9 #Found estimate of emissivity of large grains of silicon dioxide (rough approximation)- ASSUMED
mass = ((4/3)*np.pi*(radius**3))*body_density #mass of asteroid in kg
k = 1


sx = 0.0872
sy = 0
sz = -0.9962
Smag = (((sx)**2)+ ((sy)**2) + ((sz)**2))**(1/2)

def yarkovsky_effect(reb_sim):
    unit_matrix = np.array([[1, 1, 1], [1, 1, 1 ], [1, 1, 1]])
    
    v_vector = np.array([[(sp[1].vx*(1.495978707e11/31557600))], [(sp[1].vy*(1.495978707e11/31557600))], [(sp[1].vz*(1.495978707e11/31557600))]])

    r_vector = np.array([[(sp[1].x*1.495978707e11)], [(sp[1].y*1.495978707e11)], [(sp[1].z*1.495978707e11)]])
    
    hx = (r_vector[1][0]*v_vector[2][0])-(r_vector[2][0]*v_vector[1][0]) #CHECK TO MAKE SURE THIS IS RIGHT
    hy = (r_vector[2][0]*v_vector[0][0])-(r_vector[0][0]*v_vector[2][0])
    hz = (r_vector[0][0]*v_vector[1][0])-(r_vector[1][0]*v_vector[0][0])
    Hmag = (((hx)**2)+ ((hy)**2) + ((hz)**2))**(1/2)

    R1s = (1/Smag)*np.array([[0, -sz, sy],[sz, 0, -sx],[-sy, sx, 0]])
    R2s = (1/(Smag**2))*np.array([[sx**2, sx*sy, sx*sz],[sx*sy, sy**2, sy*sz],[sx*sz, sy*sz, sz**2]])
    R1h = (1/Hmag)*np.array([[0, -hz, hy],[hz, 0, -hx],[-hy, hx, 0]])
    R2h = (1/(Hmag**2))*np.array([[hx**2, hx*hy, hx*hz],[hx*hy, hy**2, hy*hz],[hx*hz, hy*hz, hz**2]])
    
    distance = (((sp[1].x**2)+(sp[1].y**2)+(sp[1].z**2))**(1/2))*1.495978707e11


    rdotv = ((r_vector[0][0]*v_vector[0][0])+(r_vector[1][0]*v_vector[1][0])+(r_vector[2][0]*v_vector[2][0]))/(c*distance)

    i_vector = (1-(rdotv))*((r_vector/(distance))-(v_vector/c)) #GONNA HAVE TO WRITE MY OWN DOT PRODUCT FUNCTION

    tanPhi = (1+(.5*(((stef_boltz*emissivity)/(np.pi**5))**(1/4))*((rotation_period/(K*C*body_density))**(1/2))*(((lsun*alph)/((distance)**2))**(3/4))))**(-1)
    tanEpsilon = (1+(.5*(((stef_boltz*emissivity)/(np.pi**5))**(1/4))*(((2*np.pi/(sp[1].n/31557600))/(K*C*body_density))**(1/2))*(((lsun*alph)/((distance)**2))**(3/4))))**(-1)
    Phi = np.arctan(tanPhi)
    Epsilon = np.arctan(tanEpsilon)

    Rys = np.dot(np.cos(Phi), unit_matrix) + np.dot(np.sin(Phi), R1s)+(np.dot((1-np.cos(Phi)), R2s))
    Ryh = np.dot(np.cos(Epsilon), unit_matrix) - (np.dot(np.sin(Epsilon),R1h))+(np.dot((1-np.cos(Epsilon)), R2h))

    yarkovsky_magnitude = (k*np.pi*(radius**2)*lsun*alph)/(4*np.pi*mass*c*((distance)**2))
    
    print(yarkovsky_magnitude)

    Yark_matrix = np.dot(Ryh, Rys)
    Direction_matrix = Yark_matrix.dot(i_vector)
    yarkovsky_acceleration = yarkovsky_magnitude*Direction_matrix
    yarkovsky_acceleration_corrected = yarkovsky_acceleration*((31557600**2)/1.495978707e11)




    sp[1].ax += yarkovsky_acceleration_corrected[0][0]				
    sp[1].ay += yarkovsky_acceleration_corrected[1][0]
    sp[1].az += yarkovsky_acceleration_corrected[2][0]

sim.additional_forces = yarkovsky_effect  # Adds force to simulation 
sim.force_is_velocity_dependent = 1 # Force depends on particle's velocity

    #Moves all particles to center of momentum frame
sim.move_to_com()
    #5.75e6
changing_a = []
changing_t = []

while sim.t < 10:                   # Max. simulation time in years
    sim.step()         # move simulation forward a time step
    changing_a.append(sp[1].a)
    print(sp[1].x)
    changing_t.append(sim.t)
    sim.integrator_synchronize()     # synchronize all changes to particles
        
plt.plot(changing_t, changing_a, label = 'Change in Semi-Major Axis')
plt.show()
print(changing_a[-1]-changing_a[0])
print((changing_a[-1]-changing_a[0])*1.495978707e11)
