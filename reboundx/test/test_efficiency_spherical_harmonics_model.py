import rebound as rb
import reboundx as rbx
import numpy as np
import time




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
    sim.add(m=1000.0, x = x_sat, y = y_sat, z = z_sat, vx = 0, vy = v, vz = 0, hash='Sat')

    return sim



def add_general_model(sim):

    rebx = rbx.Extras(sim)
    gh = rebx.load_force("gravitational_harmonics")

    J2  =  1.0826e-3
    J4  = -1.6196e-6
    J6  =  5.4068e-7
    J8  = -1.5e-7
    J10 =  5.0e-8

    N = 4
    size = (N+1)*(N+2)//2
    G = sim.G
    R_eq = 6378.137e3 #m

    # Normalized coefficients using the 4-pi normalization
    C = np.zeros(size, dtype=np.float64)
    S = np.zeros(size, dtype=np.float64)

    # The relation between the coefficients and its normalization are a /np.sqrt(2n+1), with a - for Jn.
    C[idx(2,0)] = J2        # Equal to [2,0]
    C[idx(4,0)] = J4        # Equal to [4,0]


    model = rbx.spherical_harmonics_model(
        C,
        S,
        G*m_earth,
        R_eq,
    )

    sim.particles["Earth"].params["spherical_harmonics_model"] = model
    
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

    J2 = 1.0826E-3
    J4 = -1.6196E-6
    J6  =  5.4068e-7
    J8  = -1.5e-7
    J10 =  5.0e-8
    sim.particles["Earth"].params["J2"] = J2
    sim.particles["Earth"].params["J4"] = J4
    #sim.particles["Earth"].params["J6"] = J6
    #sim.particles["Earth"].params["J8"] = J8
    #sim.particles["Earth"].params["J10"] = J10

    sim.particles["Earth"].params["R_eq"] = 6378.137e3 # m: Ecuatorial
    
    rebx.add_force(gh)

    sim._gh_force = gh
    sim._rebx = rebx





def test_consistency():

    times = np.linspace(0, 5000, 101)

    sim1 = create_sim()
    add_J2_J4(sim1)
    sim1.integrator = "ias15"
    for t in times:
        sim1.integrate(t)

    sim2 = create_sim()
    add_general_model(sim2)
    sim2.integrator = "ias15"
    for t in times:
        sim2.integrate(t)



    sat1 = sim1.particles["Sat"]
    sat2 = sim2.particles["Sat"]

    print(f' Position difference: {np.linalg.norm([sat1.x-sat2.x, sat1.y-sat2.y, sat1.z-sat2.z])}')

    assert np.allclose(
        [sat1.x, sat1.y, sat1.z],
        [sat2.x, sat2.y, sat2.z],
        rtol=1e-12,
        atol=1e-12
    )

    assert np.allclose(
    [sat1.vx, sat1.vy, sat1.vz],
    [sat2.vx, sat2.vy, sat2.vz],
    rtol=1e-12,
    atol=1e-12,
    )

    assert np.allclose(
        [sat1.ax, sat1.ay, sat1.az],
        [sat2.ax, sat2.ay, sat2.az],
        rtol=1e-13,
        atol=1e-13,
    )

    print("Everything is fine!")


    print('Accelerations:')
    print(f' For J2: {sim1.particles["Sat"].ax, sim1.particles["Sat"].ay, sim1.particles["Sat"].az}')
    print(f' For Model: {sim2.particles["Sat"].ax, sim2.particles["Sat"].ay, sim2.particles["Sat"].az}')


def test_efficiency():
    n_runs = 20
    times_J2 = []
    times_model = []
    times = np.linspace(0, 5000, 101)

    for _ in range(n_runs):
        sim = create_sim()
        add_J2_J4(sim)
        sim.integrator = "ias15"

        t0 = time.perf_counter()
        for t in times:
            sim.integrate(t)
        times_J2.append(time.perf_counter() - t0)

    for _ in range(n_runs):
        sim = create_sim()
        add_general_model(sim)
        sim.integrator = "ias15"

        t0 = time.perf_counter()
        for t in times:
            sim.integrate(t)
        times_model.append(time.perf_counter() - t0)

    mean_J2 = np.mean(times_J2)
    std_J2 = np.std(times_J2)
    mean_model = np.mean(times_model)
    std_model = np.std(times_model)

    print(f"J2/J4 legacy:    {mean_J2*1e3:.3f} +/- {std_J2*1e3:.3f} ms")
    print(f"General model:   {mean_model*1e3:.3f} +/- {std_model*1e3:.3f} ms")
    print(f"Speedup factor: {mean_model/mean_J2:.2f}x")


if __name__ == "__main__":
    test_consistency()
    test_efficiency()
