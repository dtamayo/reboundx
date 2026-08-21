import numpy as np
import time
import rebound as rb
import reboundx as rbx


def idx(n, m):
    return n * (n + 1) // 2 + m


def build_sim_with_degree(N, seed=42):
    r_earth = 6371e3
    m_earth = 5.9736e24

    sim = rb.Simulation()
    sim.units = ('m', 's', 'kg')
    sim.add(m=m_earth, r=r_earth, hash="Earth")

    r0 = 7000e3
    mu = sim.G * m_earth
    v0 = np.sqrt(mu / r0)
    sim.add(m=1000.0, x=r0, y=0, z=0, vx=0, vy=v0, vz=0, hash="Sat")

    rebx_ext = rbx.Extras(sim)
    gh = rebx_ext.load_force("gravitational_harmonics")

    size = (N + 1) * (N + 2) // 2
    rng = np.random.default_rng(seed)

    C = np.zeros(size); S = np.zeros(size)
    for n in range(2, N + 1):
        for m in range(0, n + 1):
            scale = 1e-6 / (n ** 2)
            C[idx(n, m)] = rng.normal(0, scale)
            if m > 0:
                S[idx(n, m)] = rng.normal(0, scale)

    model = rbx.spherical_harmonics_model(C, S, sim.G * m_earth, 6378.137e3)
    sim.particles["Earth"].params["spherical_harmonics_model"] = model

    sim._gh_force = gh
    sim._spherical_harmonics_model = model
    sim._C_coeffs = C
    sim._S_coeffs = S
    sim._rebx = rebx_ext

    rebx_ext.add_force(gh)
    sim.integrator = "ias15"
    return sim


def benchmark_by_degree(degrees=(2, 10, 20, 50), n_runs=10, t_max=5000, n_steps=101):
    times_arr = np.linspace(0, t_max, n_steps)
    results = {}

    for N in degrees:
        run_times = []
        for _ in range(n_runs):
            sim = build_sim_with_degree(N)
            t0 = time.perf_counter()
            for t in times_arr:
                sim.integrate(t)
            run_times.append(time.perf_counter() - t0)

        mean_t = np.mean(run_times)
        std_t = np.std(run_times)
        results[N] = (mean_t, std_t)
        print(f"N={N:3d}: {mean_t*1e3:8.3f} +/- {std_t*1e3:.3f} ms  "
              f"({mean_t/n_steps*1e6:.2f} us/step)")

    return results


if __name__ == "__main__":
    benchmark_by_degree()