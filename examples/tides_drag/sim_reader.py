import numpy as np
import matplotlib.pyplot as plt

# sun = open('Sun_Data.txt', 'r')
# sun_data = []
# for line in sun:
#     sun_data.append(line)
# sun_time = np.array([])
# sun_mass = np.array([])
# sun_radius = np.array([])

# for i in sun_data:
#     sun_time = np.append(sun_time, (i.split(('\t\t'))[0]))
#     sun_mass = np.append(sun_mass, (i.split(('\t\t'))[1]))

# sun_time = np.delete(sun_time, 0)
# sun_mass = np.delete(sun_mass, 0)
# sun_mass_corrected = np.array([])

# for i in sun_mass:
#     sun_mass_corrected = np.append(sun_mass_corrected, str(i).strip())

# sun_mass = sun_mass_corrected.astype(np.float)
# sun_time = sun_time.astype(np.float)    
# """print(sun_time)
# print(sun_mass)"""


planet = open('planet.txt', 'r')
planet_data = []
for line in planet:
    planet_data.append(line)

planet_time = np.array([])
planet_mass = np.array([])
planet_a = np.array([])
planet_e = np.array([])
planet_inc = np.array([])
planet_Omega = np.array([])
planet_omega = np.array([])
planet_f = np.array([])

for i in planet_data:
    planet_time = np.append(planet_time, (i.split(('\t\t'))[0]))
    planet_mass = np.append(planet_mass, (i.split(('\t\t'))[1]))
    planet_a = np.append(planet_a, (i.split(('\t\t'))[2]))
    planet_e = np.append(planet_e, (i.split(('\t\t'))[3]))
    planet_inc = np.append(planet_inc, (i.split(('\t\t'))[4]))
    planet_Omega = np.append(planet_Omega, (i.split(('\t\t'))[5]))
    planet_omega = np.append(planet_omega, (i.split(('\t\t'))[6]))
    planet_f = np.append(planet_f, (i.split(('\t\t'))[7]))

planet_time = np.delete(planet_time, 0)
planet_mass = np.delete(planet_mass, 0)
planet_a = np.delete(planet_a, 0)
planet_e = np.delete(planet_e, 0)
planet_inc = np.delete(planet_inc, 0)
planet_Omega = np.delete(planet_Omega, 0)
planet_omega = np.delete(planet_omega, 0)
planet_f = np.delete(planet_f, 0)

planet_f_corrected = np.array([])
for i in planet_f:
    planet_f_corrected = np.append(planet_f_corrected, str(i).strip())
    
planet_time = planet_time.astype(np.float)
planet_mass = planet_mass.astype(np.float)
planet_a = planet_a.astype(np.float)
planet_e = planet_e.astype(np.float)
planet_inc = planet_inc.astype(np.float)
planet_Omega = planet_Omega.astype(np.float)
planet_omega = planet_omega.astype(np.float)
planet_f = planet_f_corrected.astype(np.float)

# print(planet_time)
# print(planet_f)
# print(planet_a)
# print(planet_mass)

#Plots the semi-major axis of the planet and the radius of the Sun over the course of the simulation
fig, ax = plt.subplots()
# plt.ylim(0, 1.5)
plt.title(r'$R_{\odot}$' ' & ' r'$a_{\oplus}$' ' vs. Simulation Time')
plt.xlabel('sim.t [yr]')
plt.ylabel('[AU]')
ax.plot(planet_time, planet_a, label = '$a_{\oplus}$')
# ax.plot(np.linspace(START_AGE, START_AGE + sim.t ,len(changing_a)), csr(np.linspace(START_AGE, START_AGE + sim.t, len(changing_a))), label = "Sun's surface radius")
ax.legend(loc='best', ncol=1)
plt.show()

# #Plots planet's eccentricity vs. sim time
# fig, ax2 = plt.subplots()
# # plt.ylim(0, 1.5)
# plt.title(r'$e_{\oplus}$' ' vs. Simulation Time')
# plt.xlabel('sim.t [yr]')
# plt.ylabel(r'$e$')
# ax2.plot(planet_time, planet_e, label = '$e_{\oplus}$')
# ax2.legend(loc='best', ncol=1)
# plt.show()