import rebound as rb
import reboundx as rbx
import numpy as np




r_earth = 6371e3 # m
m_earth = 5.9736*10**24 # kg

def idx(n, m):
    return n*(n+1)//2 + m


def create_sim():
    sim = rb.Simulation()
    sim.units = ('m', 's', 'kg') 

    sim.add(m=m_earth, r=r_earth, hash=str("Earth"))
    
    lon = np.deg2rad(45)
    lat = np.deg2rad(45)

    r = 7000e3
    x_sat = r*np.cos(lat)*np.sin(lon)
    y_sat = r*np.cos(lat)*np.cos(lon)
    z_sat = r*np.sin(lat)
    mu = sim.G * m_earth
    v = np.sqrt(mu/r)
    sim.add(m=10.0, x = x_sat, y = y_sat, z = z_sat, vx = 0, vy = v, vz = 0, hash='Sat')

    return sim



def add_general_model(sim):

    rebx = rbx.Extras(sim)
    gh = rebx.load_force("gravitational_harmonics")

    J2  =  1.0826e-3
    J4  = -1.6196e-6
    J6  =  5.4068e-7
    J8  = -1.5e-7
    J10 =  5.0e-8

    N = 50
    size = (N+1)*(N+2)//2
    G = sim.G
    R_eq = 6378.137e3 #m

    # Normalized coefficients using the 4-pi normalization
    C = np.zeros(size, dtype=np.float64)
    S = np.zeros(size, dtype=np.float64)

    # The relation between the coefficients and its normalization are a /np.sqrt(2n+1), with a - for Jn.
    C[idx(2,0)] = -J2/np.sqrt(5) # Equal to [2,0]
    C[idx(4,0)] = -J4/3          # Equal to [4,0]
    C[idx(6,0)] = -J6/np.sqrt(13)
    C[idx(8,0)] = -J8/np.sqrt(17)
    C[idx(10,0)] = -J10/np.sqrt(21)



    model = rbx.spherical_harmonics_model(
        C,
        S,
        G*m_earth,
        R_eq,
    )

    sim.particles["Earth"].params["spherical_harmonics_model"] = model
    sim.particles["Earth"].params["Omega"] = [0.0, 0.0, 7.292115e-5]
    
    # Given the model made in a function, we need to save different elements as long as the simulation lasts
    sim._gh_force = gh
    sim._spherical_harmonics_model = model
    sim._C_coeffs = C   # if the model only saves the pointer to buffer numpy
    sim._S_coeffs = S
    sim._rebx = rebx

    rebx.add_force(gh)



def add_J2_J4(sim):

    rebx = rbx.Extras(sim)
    gh = rebx.load_force("gravitational_harmonics")

    J2  =  1.0826e-3
    J4  = -1.6196e-6
    J6  =  5.4068e-7
    J8  = -1.5e-7
    J10 =  5.0e-8

    sim.particles["Earth"].params["J2"] = J2
    sim.particles["Earth"].params["J4"] = J4
    sim.particles["Earth"].params["J6"] = J6
    sim.particles["Earth"].params["J8"] = J8
    sim.particles["Earth"].params["J10"] = J10
    sim.particles["Earth"].params["Omega"] = [0.0, 0.0, 7.292115e-5]    
    sim.particles["Earth"].params["R_eq"] = 6378.137e3 # m: Ecuatorial radius
    
    rebx.add_force(gh)





def test_compatibility():

    times = np.linspace(0, 15000, 1001)

    sim1 = create_sim()
    sim2 = create_sim()

    add_J2_J4(sim1)
    add_general_model(sim2)

    sim1.integrator = "ias15"
    sim2.integrator = "ias15"  

    print('Integrating...')

    for time_integrate in times: 
        

        sim1.integrate(time_integrate)
        sim2.integrate(time_integrate)

        sat1 = sim1.particles["Sat"]
        sat2 = sim2.particles["Sat"]


        rtol = 1e-6
        atol = 1e-6

        assert np.allclose(
            [sat1.x, sat1.y, sat1.z],
            [sat2.x, sat2.y, sat2.z],
            rtol=rtol,
            atol=atol
        )

        assert np.allclose(
        [sat1.vx, sat1.vy, sat1.vz],
        [sat2.vx, sat2.vy, sat2.vz],
        rtol=rtol,
        atol=atol,
        )

        assert np.allclose(
            [sat1.ax, sat1.ay, sat1.az],
            [sat2.ax, sat2.ay, sat2.az],
            rtol=rtol,
            atol=atol,
        )


    print('Accelerations:')
    print(f' For J2: {sim1.particles["Sat"].ax, sim1.particles["Sat"].ay, sim1.particles["Sat"].az}')
    print(f' For Model: {sim2.particles["Sat"].ax, sim2.particles["Sat"].ay, sim2.particles["Sat"].az}')

    print(f'Difference in precision xyz: {sim1.particles["Sat"].x - sim2.particles["Sat"].x, sim1.particles["Sat"].y - sim2.particles["Sat"].y, sim1.particles["Sat"].z - sim2.particles["Sat"].z}')


test_compatibility()
