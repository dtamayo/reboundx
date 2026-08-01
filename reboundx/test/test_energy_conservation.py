import numpy as np
import rebound as rb
import reboundx as rbx
from test_utils import spherical_harmonics_value


def idx(n, m):
    return n * (n + 1) // 2 + m


def total_energy(sim, model):

    E_newtonian = sim.energy()

    earth = sim.particles["Earth"]
    sat = sim.particles["Sat"]

    dx, dy, dz = sat.x - earth.x, sat.y - earth.y, sat.z - earth.z
    r = np.sqrt(dx*dx + dy*dy + dz*dz)
    phi = np.arcsin(dz / r)
    lam = np.arctan2(dy, dx)


    U_pert = spherical_harmonics_value(model, phi, r, lam)
    E_pert = -sat.m * U_pert 
    return E_newtonian + E_pert


def test_energy_conservation_long_orbit():
    r_earth = 6371e3
    m_earth = 5.9736e24

    sim = rb.Simulation()
    sim.units = ('m', 's', 'kg')
    sim.add(m=m_earth, r=r_earth, hash="Earth")

    r0 = 7000e3
    mu = sim.G * m_earth
    v0 = np.sqrt(mu / r0)
    sim.add(m=1000.0, x=r0, y=0, z=0, vx=0, vy=v0, vz=v0*0.3, hash="Sat")  

    rebx_ext = rbx.Extras(sim)
    gh = rebx_ext.load_force("gravitational_harmonics")

    N = 15
    size = (N + 1) * (N + 2) // 2
    C = np.zeros(size); S = np.zeros(size)
    J2, J4 = 1.0826e-3, -1.6196e-6
    C[idx(2, 0)] = J2 
    C[idx(4, 0)] = J4 

    model = rbx.spherical_harmonics_model(C, S, sim.G * m_earth, 6378.137e3)
    sim.particles["Earth"].params["spherical_harmonics_model"] = model

    sim._gh_force = gh
    sim._spherical_harmonics_model = model
    sim._C_coeffs = C
    sim._S_coeffs = S
    sim._rebx = rebx_ext

    rebx_ext.add_force(gh)
    sim.integrator = "ias15"

    E0 = total_energy(sim, model)

    n_orbits = 50
    T = 2 * np.pi * np.sqrt(r0**3 / mu)
    times = np.linspace(0, n_orbits * T, 500)

    max_rel_drift = 0.0
    for t in times:
        sim.integrate(t)
        E = total_energy(sim, model)
        rel_drift = abs((E - E0) / E0)
        max_rel_drift = max(max_rel_drift, rel_drift)

    print(f"Max relative energy drift over {n_orbits} orbits: {max_rel_drift:.3e} Jules")
    assert max_rel_drift < 1e-8, f"Energy drift too large: {max_rel_drift} Jules"


if __name__ == "__main__":
    test_energy_conservation_long_orbit()