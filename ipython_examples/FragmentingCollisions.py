import rebound
import reboundx
import matplotlib.pyplot as plt
import numpy as np

sim = rebound.Simulation()
rebx = reboundx.Extras(sim)
# Load "fragmenting_collisions" as the collision resolve module
collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
rebx.add_collision_resolve(collision_resolve)
sim.integrator = "mercurius"
sim.collision = "direct"
#sim.units = ('yr', 'AU', 'Msun')
sim.G = 39.476926421373
sim.dt = 6.0/365.0
sim.rand_seed = 1

sim.add(m = 1, r = 0.00465)

n_pl = 30 # Number of planetesimals

# planetesimal mass range
lunar_mass  = 3.8e-8    # solar masses
earth_mass  = 3e-6      # solar masses
mass_min    = 0.5 * lunar_mass
mass_max    = 0.1 * earth_mass

# Set up minimum fragment mass
collision_resolve.params["fc_min_frag_mass"] = mass_min
collision_resolve.params["fc_particle_list_file"] = "family_tree_python.csv"

rho = 5.05e6   # units of solMass/AU^3, equal to ~3 g/cm^3

# Add 30 planetary embryos
for i in range(n_pl):
    # Orbital parameters
    a = np.random.uniform(0.1, 0.5) # semi-major axis in AU
    e = np.random.uniform(0.0, 0.01) # eccentricity
    inc = np.random.uniform(0.0, np.pi/180) # inclination in radians
    omega = np.random.uniform(0.0, 2.*np.pi) # argument of periapsis
    Omega = np.random.uniform(0.0, 2.*np.pi) # longitude of ascending node
    f = np.random.uniform(0.0, 2.*np.pi) # mean anomaly
    m = np.random.uniform(mass_min, mass_max) # particle mass in solar masses
    r = ((3*m)/(4 * np.pi * rho)) ** (1/3) * 10 # particle radius in AU, times 10 to ease collisions
    sim.add(m=m, a=a, e=e, inc=inc, omega=omega, Omega=Omega, f=f, r=r)

sim.move_to_com()

integration_time = 1e3

sim.integrate(integration_time)