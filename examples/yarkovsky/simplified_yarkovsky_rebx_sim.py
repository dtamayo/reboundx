import rebound
import reboundx
import numpy as np

#Simulation begins here
sim = rebound.Simulation()

sp = sim.particles #simplifies way to access particles parameters 

sim.units = ('yr', 'AU', 'Msun') #changes simulation and G to units of solar masses, years, and AU  
sim.integrator = "whfast" #integrator for sim
sim.dt = .1 #timestep for sim


#Adds Sun and outer planets to simulation
sim.add(m=1) 
sim.add(a=.5, f=0, Omega=0, omega=0, e=0, inc=0, m=0)

#Moves all particles to center of momentum frame
sim.move_to_com()

#Gives orbital information before the simulation begins
print("\n***INITIAL ORBITS:***")
for orbit in sim.calculate_orbits():
    print(orbit)
    
rebx = reboundx.Extras(sim)
yark = rebx.load_force("max_yarkovsky")

sp[1].params["body_density"] = 3000
sp[1].params["direction_flag"] = 1
yark.params["lstar"] = 3.828e26
sp[1].r = 1000 

rebx.add_force(yark)

print(sp[1].params['body_density'])

tmax=50000 # in yrs


a = .5 #starting semi-major axis for the asteroid
            
sim.integrate(tmax) #integrates system for tmax years       
 
final_a = sp[1].a #semi-major axis of asteroid after the sim    
                      
print("CHANGE IN SEMI-MAJOR AXIS", final_a-a) #prints difference between the intitial and final semi-major axes of asteroid
