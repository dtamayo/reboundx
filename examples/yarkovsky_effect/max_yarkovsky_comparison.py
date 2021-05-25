import rebound
import reboundx
import numpy as np
import matplotlib.pyplot as plt

#Simulation begins here
sim = rebound.Simulation()

sp = sim.particles #simplifies way to access particles parameters 

sim.units = ('yr', 'AU', 'Msun') #changes simulation and G to units of solar masses, years, and AU  
sim.integrator = "mercurius" #integrator for sim
sim.dt = .5 #timestep for sim

sim.add(m=1) #Adds Sun 
sim.add(a=.5, f=0, Omega=0, omega=0, e=0, inc=0, m=0) #Adds test particle 

#Moves all particles to center of momentum frame
sim.move_to_com()

#Gives orbital information before the simulation begins
print("\n***INITIAL ORBITS:***")
for orbit in sim.calculate_orbits():
    print(orbit)
    
rebx = reboundx.Extras(sim)
yark = rebx.load_force("max_yarkovsky")

au_conv = 1.495978707e11
msun_conv = 1.9885e30
yr_conv = 31557600.0

#Converts units of parameters from m/kg/sec to AU/Msun/yr
my_density = (3000*au_conv*au_conv*au_conv)/msun_conv
my_c = (2.998e8*yr_conv)/au_conv
my_lstar = (100000*3.828e26*yr_conv*yr_conv*yr_conv)/(msun_conv*au_conv*au_conv)
radius = 1000/au_conv


yark.params["my_lstar"] = my_lstar
yark.params["my_c"] = my_c
sp[1].params["my_body_density"] = my_density
sp[1].params["direction_flag"] = 1
sp[1].r = radius

rebx.add_force(yark)

tmax=50000 # in yrs

a_start = .5 #starting semi-major axis for the asteroid
        
changing_a = []
    
while sim.t < tmax: #determines how many years the simulation goes for
        sim.step() #moves simulation forward a time step
        sim.integrator_synchronize() #synchronizes all changes in the particles during the simulation
        if (sim.t % 1000) < 0.05: # every 1000 yr
            changing_a.append(sp[1].a) #adds semi-major axis to list every thousand years      
 
a_final = sp[1].a #semi-major axis of asteroid after the sim    
                      
print("CHANGE IN SEMI-MAJOR AXIS", a_final-a_start, "AU") #prints difference between the intitial and final semi-major axes of asteroid

print(sp[1].z)


##########################################




#Simulation begins here
sim = rebound.Simulation()

sp = sim.particles #simplifies way to access particles parameters 

sim.units = ('yr', 'AU', 'Msun') #changes simulation and G to units of solar masses, years, and AU  
sim.integrator = "whfast" #integrator for sim
sim.dt = .05 #timestep for sim

sim.add(m=1) #Adds Sun 
sim.add(a=.5, f=0, Omega=0, omega=0, e=0, inc=0, m=0) #Adds test particle 

#Moves all particles to center of momentum frame
sim.move_to_com()

#Gives orbital information before the simulation begins
print("\n***INITIAL ORBITS:***")
for orbit in sim.calculate_orbits():
    print(orbit)
    
au_conv = 1.495978707e11
msun_conv = 1.9885e30
yr_conv = 31557600.0

#Converts units of parameters from m/kg/sec to AU/Msun/yr
my_density = (3000*au_conv*au_conv*au_conv)/msun_conv
my_c = (2.998e8*yr_conv)/au_conv
my_lstar = (100000*3.828e26*yr_conv*yr_conv*yr_conv)/(msun_conv*au_conv*au_conv)
radius = 1000/au_conv

def yark_approx(reb_sim):
    
    distance = ((sp[1].x**2)+(sp[1].y**2)+(sp[1].z**2))**(1/2)
    
    speed = ((sp[1].vx**2)+(sp[1].vy**2)+(sp[1].vz**2))**(1/2)
    
    yarkovsky_magnitude = (my_lstar)/(4*((4*np.pi*radius*(my_density))/3)*(my_c)*distance*distance) #FIX ME
    
    sp[1].ax += yarkovsky_magnitude*(sp[1].vx/speed)
    sp[1].ay += yarkovsky_magnitude*(sp[1].vy/speed)
    sp[1].az += yarkovsky_magnitude*(sp[1].vz/speed)
    
sim.additional_forces = yark_approx  # adds force to simulation 
sim.force_is_velocity_dependent = 1

tmax=50000 # in yrs

a_start = .5 #starting semi-major axis for the asteroid

changing_a_approx = []    
                    
while sim.t < tmax: #determines how many years the simulation goes for
#        start = time.perf_counter() #beginning of time step timer
        sim.step() #moves simulation forward a time step
        sim.integrator_synchronize() #synchronizes all changes in the particles during the simulation
        if (sim.t % 1000) < 0.05: # every 1000 yr
            changing_a_approx.append(sp[1].a) #adds semi-major axis to list every thousand years #integrates system for tmax years       
 
a_final = sp[1].a #semi-major axis of asteroid after the sim    
                      
print("CHANGE IN SEMI-MAJOR AXIS", a_final-a_start, "AU") #prints difference between the intitial and final semi-major axes of asteroid

#Uses anaylytical equation to plot the change in semi-major axis
G = 4*np.pi**2
mass = (4/3)*np.pi*radius**3*my_density


time = np.linspace(0, sim.t, len(changing_a))

a = []

for i in time:
    dist = ((3/((G*1)**(1/2)))*((radius**2*my_lstar*i)/(4*mass*my_c))+a_start**(3/2))**(2/3)
    a.append(dist)
    
a_2 = []

for i in time:
    dist = ((3*radius**2*my_lstar*i)/(32*2*np.pi*mass*my_c)+a_start**(3/2))**(2/3)
    a_2.append(dist)

"""fig1, ax1 = plt.subplots()
#plt.axis([0, sim.t, 1.126016, 1.126015])
plt.title('Max Yarkovsky Comparison')
plt.xlabel('Simulation Time [yr]')
plt.ylabel('[AU]')
#ax1.set_yscale('log')
ax1.plot(np.linspace(0, sim.t, len(changing_a)),
         changing_a, label = 'Max Yarkovsky Rebx Sim')
ax1.plot(np.linspace(0, sim.t, len(changing_a_approx)),
         changing_a_approx, label = 'Yarkovsky Approximation Sim')
ax1.plot(time, a, label = 'My Yarkovsky Equation')
ax1.plot(time, a_2, label = 'Veras 2015 Yarkovsky Equation')
plt.legend()

plt.show()"""

print(a[-1]/a_2[-1])
print(changing_a_approx[-1]/changing_a[-1])

fig, ax1 = plt.subplots()
ax1.plot(np.linspace(0, sim.t, len(changing_a)),
         changing_a, '-', color = 'red', label = 'Rebound Sim w/ Yarkovsky Rebx Effect')
ax1.plot(time, a_2, ':', color = 'grey', label = 'Veras et. al (2015) Equation')
ax1.set_xlabel('Simulation Time (yr)')
ax1.set_ylabel('Asteroid Semi-major Axis (AU)')
ax1.set_yscale('log')
ax1.legend()

plt.show()