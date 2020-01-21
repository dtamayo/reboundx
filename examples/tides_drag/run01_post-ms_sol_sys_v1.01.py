# =============================================================================
# Post-MS Solar System Orbital Evolution (REBOUND)
# - TidalForce() Optimization Factor: 5000x
# - Coefficient spline for stellar-dependent variables (smallGamma)
# - Checkpoint / 0.5 Myr
# - OMEGAsun = 0 (Sun's surface ang. vel.)
# - Includes all eight planets
# - ~300 Myr total sim time
# - Planet's mass goes to zero when engulfed
# - Data saved to text files
#
# Authors: Noah Ferich, Stanley A. Baronett
# Last Updated: 8/28/2019
# Version: 1.01
# =============================================================================

import mesa_reader as mr
import time
import rebound
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import os

# =============================================================================
# Import MESA history.data using PyMesaReader
# =============================================================================
# Make MesaData objects from a history files
h = mr.MesaData('LOGS/history_to_end_agb.data')
h1 = mr.MesaData('LOGS/history_to_wd.data')

# Input sequence of data files into a, m, r, l and omega arrays

# SOLAR AGE DATA (measured in years)
# (First real age entry at approx. 1.2275e+10 yrs;
# initial 3 points used for smooth transition of 10 millions years)
a = np.array([1.226500e+10, 1.226525e+10, 1.226550e+10])
a1 = h.star_age
a2 = h1.star_age

# SOLAR MASS DATA (measured in Msun units)
# (First 3 data points used for smooth transition from 1 Msun)
m = np.array([1.0, 1.0, 1.0]) 
m1 = h.star_mass
m2 = h1.star_mass

# SOLAR RADIUS DATA (measured in RsunAUAU units)
# (First 3 data points used for smooth transition from 1 RsunAUAU)
r = np.array([1.0, 1.0, 1.0]) 
r1 = h.radius
r2 = h1.radius

# SOLAR LUMINOSITY DATA (measured in Lsun units)
# (First 3 data points used for smooth transition from 1 Lsun)
l = np.array([1.0, 1.0, 1.0]) 
l1 = h.luminosity
l2 = h1.luminosity

# Concatenate data into single arrays
a = np.append(a, a1)
m = np.append(m, m1)
r = np.append(r, r1)
l = np.append(l, l1)
age = np.append(a, a2)
Msun = np.append(m, m2)
Rsun = np.append(r, r2)
Lsun = np.append(l, l2)

# Convert Rsun values into AU (REBOUND default units)
RsunAU = np.zeros(Rsun.size)
for i, rad in enumerate(Rsun):
    RsunAU[i] = rad * 0.00465047 # approx. 215 Rsun to 1 AU

# Convert LSun into watts (MKS units)
Lsun_Watts = np.zeros(Lsun.size)
for i, l  in enumerate(Lsun):
    Lsun_Watts[i] = l * 3.828e26 # IAU Resolution B3 conversion

#Convert LSun_Watts into REBOUND AU, Msun, and yrs
Lsun_Sim_Units = np.zeros(Lsun.size)
for i, w in enumerate(Lsun_Watts):
    Lsun_Sim_Units[i] = (w * ((6.7e-12)**2) * (5e-31)) / ((3.2e-8)**3) 

csm = CubicSpline(age, Msun)                 # mass function  
csr = CubicSpline(age, RsunAU)               # radius function
csl = CubicSpline(age, Lsun_Sim_Units)       # luminosity function
#csOmega = CubicSpline(age, OMEGAsunYR)       # angular velocity function

# =============================================================================
# RECREATE SCHRODER DATA: spline functions with natural boundary condition type
# =============================================================================
"""schroder_age = np.array([0.0, 2.6e6, 3.55e6, 3.925e6])
schroder_mass_Msun = np.array([1.0, 0.9, 0.8, 0.7])
schroder_radius_AU = np.array([0.4, 0.7, 0.966, 1.16])
cs_natural_r = CubicSpline(schroder_age, schroder_radius_AU, bc_type='natural')
cs_natural_m = CubicSpline(schroder_age, schroder_mass_Msun, bc_type='natural')

cs_natural_l = CubicSpline(age, Lsun_Sim_Units, bc_type='natural')"""

# =============================================================================
# Small gamma coefficient precalculation:
# All stellar-dependent variables are calculated for the entire life of the
# star (from MESA data), stored in an array and splined before the simulation
# starts and TidalForce() is applied.
# =============================================================================
#small_gamma_time = np.linspace(0, age, age.all())

# Convective friction time
t_f = np.zeros(age.size)
for i in range(len(t_f)):
    t_f[i] = (csm(age[i]) * (csr(age[i])**2) / csl(age[i]))**(1/3)
    
# Coefficient lambda_2 depends on the properties of the convective envelope
# From Schroder & Smith (2008), "For a fully convective envelope (Zahn 1989,
# Eq. 15), with a tidal period of ~ O(1y), comparable to 2t_f, we may use
# lambda_2 ~ 0.019 * alpha^(4/3) ~ 0.038 (with a convection parameter of our
# tip-RGB solar model of alpha ~ 1.7)"
lambda_2 = 0.038
# Precalculate small gamma values, i.e., all stellar-dependent variables into
# one coefficient used later in Gamma (tidal torque).
smallGamma = np.zeros(age.size)
for i in range(len(smallGamma)):
    smallGamma[i] = 6 * (lambda_2 / t_f[i]) * (1/csm(age[i]))**2 * csm(age[i]) * (csr(age[i])**8)

cs_smallGamma = CubicSpline(age, smallGamma) # coefficients for TidalTorque()

# =============================================================================
# Restart existing REBOUND Simulation
# =============================================================================
#Creates new simulation using checkpoint information
"""sim = rebound.Simulation(os.path.join(chkptpath, 'Simulation Checkpoint-##########years'))
sim.status()"""

# =============================================================================
# Creating REBOUND Simulation
# =============================================================================
sim = rebound.Simulation()

# Time step for the simulation
sim.dt = 0.0120423 # 1/20 Mercury's orbital period
#sim.dt = .05 # 1/20 Earth's orbital period
sim.integrator = "mercurius" # changes which integrator the simulation uses
sp = sim.particles

#Changes simulation and G to units of solar masses, years, and AU
sim.units = ('yr', 'AU', 'Msun') 
#START_AGE = 12389500000   # Schroder test beginning
START_AGE = age[0]         # Where simulation begins
TIP_RGB_AGE = 12393500000  # in years

#Adds Sun and outer planets to simulation
#Sun starts with the mass at START_AGE 
sim.add(m=1, x=0, y=0, z=0, hash = 'sun') 
sim.add(m=1.6601141530543488e-07, a=0.38709830633484105, e=0.20563990000064852, inc=0.12223965867367434, Omega=0.8431133962407107, omega=0.5093185547138174, f=-0.9908892789276185, hash='mercury')
sim.add(m=2.4478382877847715e-06, a=0.7385168404826685, e=0.01791368562875911, inc=0.058930414887247075, Omega=1.340369994329905, omega=1.2600802146822982, f=-0.24227419974291098, hash='venus' )
sim.add(m=3.040432648022642e-06, a=0.9857355919128629, e=0.022089569577719258, inc=2.9693034327917464e-05, Omega=3.1174816544358386, omega=-0.9862996632438596, f=3.416281877452467, hash='earth')
sim.add(m=3.2271560375549977e-07, a=1.5348433917046136, e=0.08919918408697937, inc=0.03211689542906033, Omega=0.8633257005591881, omega=-1.3024394497901224, f=3.045949491277054, hash='mars')
sim.add(m=0.0009547919152112404, a=5.177672552407451, e=0.05042022658304223, inc=0.02274728959386692, Omega=1.7541044834434427, omega=-1.4438820758105815, f=4.298882445213031, hash='jupiter')
sim.add(m=0.0002858856727222417, a=9.501512380472281, e=0.05814993650690009, inc=0.04343882978676329, Omega=1.9827226495964432, omega=-0.3513794685846589, f=3.3937790865431725, hash='saturn')
sim.add(m=4.36624373583127e-05, a=19.23955552089647, e=0.04662499096808216, inc=0.013500700809991764, Omega=1.2904337179090757, omega=1.635605857481013, f=-2.3411557180758407, hash='uranus')
sim.add(m=5.151383772628674e-05, a=30.024914413466103, e=0.010211207912002765, inc=0.030889408668620208, Omega=2.300001316769195, omega=-1.2692861280628678, f=-1.259521330581421, hash='neptune')

#Planets to have tidal force applied to them
Tidal_Planets = [1, 2, 3]
tidal_planet_index = 0
current_tidal_planet = Tidal_Planets[tidal_planet_index]

#Moves all particles to center of momentum frame
sim.move_to_com()

# =============================================================================
# Defining Helper Functions
# =============================================================================
# Function for finding the distance to a planet's pericenter
def Pericenter(a, e):
    return a * (1-e)

# Function for finding Distance to a planet 
def PlanetDistance(i):
        return (sp[i].x**2+sp[i].y**2+sp[i].z**2)**(1/2)      

# Torque on planet from tidal bulge of Sun                                             
def TidalTorque():  
    return (sp[current_tidal_planet].m)**2*(cs_smallGamma(START_AGE + sim.t)/PlanetDistance(current_tidal_planet)**6)*(0 - (2*np.pi/sp[current_tidal_planet].P))

# Main function added to simulation for the tidal force
stop_flag = False 
def TidalPlanet(): # Determines what planet will be affected by the tidal force
    global current_tidal_planet
    global tidal_planet_index
    global stop_flag
    if sp[current_tidal_planet].m == 0:
            tidal_planet_index += 1
            if tidal_planet_index == len(Tidal_Planets):
                stop_flag = True
            if stop_flag == False:
                current_tidal_planet = Tidal_Planets[tidal_planet_index]

stop_flag_2 = False        
closest_planet = 1
def FirstPlanet(): #Determines the closest planet
    global stop_flag_2
    global closest_planet
    if stop_flag_2 == False:
        if sp[closest_planet].m == 0:
            closest_planet += 1 
        if closest_planet == len(sp):
            stop_flag_2 = True
FirstPlanet() #Initializes before simulation
        
OPTIMIZATION_FACTOR = 5000
Step_Counter = 1
def TidalForce(reb_sim):
    global Step_Counter
    global current_tidal_planet
    if stop_flag == False:        
            TidalPlanet()
            if  Step_Counter % OPTIMIZATION_FACTOR == 0:
                # Calculate & store values needed for accel. change
                torque = OPTIMIZATION_FACTOR * TidalTorque()
                p_distance = PlanetDistance(current_tidal_planet)
                p_mass = sp[current_tidal_planet].m
                p_total_vel = (sp[current_tidal_planet].vx**2+sp[current_tidal_planet].vy**2+sp[current_tidal_planet].vz**2)**(1/2)
                # Adjust accel. vectors of planet
                sp[current_tidal_planet].ax += (torque/(p_mass*p_distance))*(sp[current_tidal_planet].vx/p_total_vel)				
                sp[current_tidal_planet].ay += (torque/(p_mass*p_distance))*(sp[current_tidal_planet].vy/p_total_vel)
                sp[current_tidal_planet].az += (torque/(p_mass*p_distance))*(sp[current_tidal_planet].vz/p_total_vel)				
                Step_Counter = 0 # Reset counter
            Step_Counter += 1


sim.additional_forces = TidalForce  # Adds force to simulation 
sim.force_is_velocity_dependent = 1 # Force depends on particle's velocity

# list that shows how the semi-major axis of planet has changed after simulation 
changing_a = []
engulfed_date = 0.0 # To store the Sun's age in years should Earth be engulfed

#Creates new files that will store simulation data
for i in range(len(sp)):
    if i == 0:
        sun_file = open("Data/Sun_Data.txt", "w+")
        sun_file.write("Time(yrs)\t\tMass(Msun)\t\tRadius(AU)\t\tLuminosity(Msun-yrs-AU)\n")
        sun_file.close()
        
    if i > 0:
        planet_file = open("Data/Planet_" + str(i) + "_Data.txt", "w+")
        planet_file.write("Time(yrs)\t\tMass(Msun)\t\tSemi-major Axis(AU)\t\tEccentricity\t\tInclination(Radians)\t\tLongitude_of_Ascending_Node(Radians)\t\tArgument_of_Periapsis(Radians))\t\tTrue_Anomaly(Radians)\n")
        planet_file.close()
# =============================================================================
# Starting Simulation
# =============================================================================
# Gives orbital information before the simulation begins
print("\n***INITIAL ORBITS:***")
for orbit in sim.calculate_orbits():
    print(orbit)
print("\n\n***SIMULATION TRACKING:***")
timer_start = time.perf_counter() #Start of timer
# Set checkpoint path
chkptpath = 'Checkpoints/'

t_prev_sim = 0           # Checks data_interval has passed
t_prev_check = 0         # Checks chkpt_interval has passed
data_interval = 1000     # in yr
chkpt_interval = 500000  # in yr

# Main simulation loop
while sim.t < (age[-1] - START_AGE): # Max. simulation time in years
    sim.step()                       # move simulation forward a time step
    sim.integrator_synchronize()     # synchronize all changes to particles
    sp[0].m = csm(sim.t + START_AGE) # adjust mass of the Sun
    if (sim.t - t_prev_sim) >= data_interval or sim.t == sim.dt:        # every 1000 yr
        t_prev_sim = sim.t
        if stop_flag_2 == False:
            changing_a.append(sp[closest_planet].a)   # Add semi-major axis to list    
            print("@ {0:1.0f}".format(sim.t), "yr: Closest Planet-a:{0:1.12f}".format(sp[closest_planet].a),
                "Rsun:{0:1.12f}".format(csr(sim.t + age[0])),
                "Msun:{0:1.12f}".format(csm(sim.t + age[0])))
        # Check planet's pericenter is within solar surface radius 
            if Pericenter(sp[closest_planet].a, sp[closest_planet].e) <= csr(sim.t + START_AGE):
                sp[0].m += sp[closest_planet].m
                sp[closest_planet].m = 0
                engulfed_date = sim.t + START_AGE
                print("!!! CLOSEST PLANET ENGULFED AT", engulfed_date, "YEARS !!!")
                FirstPlanet()
        for i in range(len(sp)): #Prints data to Sun and planets
            if i == 0:
                sun_file = open("Data/Sun_Data.txt", "a+")
                sun_file.write(str(sim.t) + "\t\t" + str(sp[0].m) + "\t\t" + str(csr(sim.t + START_AGE)) + "\t\t" + str(csl(sim.t + START_AGE)) + "\n")
                sun_file.close()
            if i > 0:
                planet_file = open("Data/Planet_" + str(i) + "_Data.txt", "a+") 
                if sp[i].m > 0: 
                    planet_file.write(str(sim.t) + "\t\t" + str(sp[i].m) + "\t\t" + str(sp[i].a) + "\t\t" + str(sp[i].e) + "\t\t" + str(sp[i].inc) + "\t\t" + str(sp[i].Omega) + "\t\t" + str(sp[i].omega) + "\t\t" + str(sp[i].f) + "\n")
                else:
                    planet_file.write(str(sim.t) + "\t\t" + str(0) + "\t\t" + str(0) + "\t\t" + str(0) + "\t\t" + str(0) + "\t\t" + str(0) + "\t\t" + str(0) + "\t\t" + str(0) + "\n")
                planet_file.close()    
    # Save checkpoint every 0.5 Myr of simulation
    if (sim.t - t_prev_check) >= chkpt_interval:
        t_prev_check = sim.t
        sim.save(os.path.join(chkptpath, "checkpoint_"+ str(int(sim.t)) + "_yr"))
        print("SIMULATION SAVED")        

# =============================================================================
# Post Simulation
# =============================================================================
# End of simulation timer and prints time in seconds
timer_end = time.perf_counter() 
runtime = timer_end - timer_start
runtime_min = runtime / 60
runtime_hr = runtime_min / 60

print("\n***SIMULATION SUMMARY:***")
print("Solar system final age:", (sim.t + START_AGE)/1e9, "Gyr")
print("Total simulation time:", sim.t/1e6, "Myr")
print("Total runtime:", runtime, "sec")
print("              ", runtime_min, "min")
print("              ", runtime_hr, "hr")
print("Sun's final mass:", sim.particles[0].m, "Msun") 
print("Last engulfed at:", engulfed_date, "yr")

# Print information on the orbits after the simulation
print("\n***FINAL ORBITS:***")
for orbit in sim.calculate_orbits():
    print(orbit)

#Plots the orbits of the planets
#fig2 = rebound.OrbitPlot(sim, trails=True, unitlabel ="[AU]")

# Plot Rsun & closest planet's semi-major axis against sim. time
if engulfed_date == 0: # If never engulfed
    engulfed_date = TIP_RGB_AGE # Make pre-tip-RGB engulfed be 0
engulfed_before_tip_RGB = TIP_RGB_AGE - engulfed_date

textstr = '\n'.join((
    r'$\mathrm{Runtime}:%.0f$ min / $%0.1f$ hr' % 
    (runtime_min, runtime_hr),
	'$M_{\odot final}=%.3f M_{\odot}$' %
	sim.particles[0].m,
	'Pre-tip-RGB engulf: $%.0f$ yr' %
	engulfed_before_tip_RGB))
fig, ax = plt.subplots(figsize=(8, 4.5))
plt.axis([START_AGE, START_AGE + sim.t, 0, 1.4])
plt.title(r'$R_{\odot}$' ' & ' r'$r_{CP}$' ' vs. Simulation Time')
plt.xlabel('Simulation Time [yr]')
plt.ylabel('[AU]')
ax.plot(np.linspace(START_AGE, START_AGE + sim.t, len(changing_a)),
         changing_a, label = r'$r_{CP}$')
ax.plot(np.linspace(START_AGE, START_AGE + sim.t ,len(changing_a)),
         csr(np.linspace(START_AGE, START_AGE + sim.t, len(changing_a))),
         label = r'$R_{\odot}$')
ax.legend(loc='best', ncol=1)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.25)
ax.text(0.02, 0.22, textstr, transform=ax.transAxes, fontsize=10,
         verticalalignment='top', bbox=props)
plt.savefig('Graphs/Engulf_Graph.png', dpi=400)
#plt.show()