import numpy as np
import time
import rebound 
import reboundx 
import random

AU = 1.495978707e11
Msun = 1.98847e30
Mearth = 5.9722e24
Rsun = 6.957e8
Rearth = 6.371e6
c_light = 299792458.0
 
def build_bench_sim(N, shadow_model, n_occulters=1):
    sim = rebound.Simulation()
    sim.units = ('m', 's', 'kg')
    sim.add(m=Msun, x=0, y=0, z=0, hash="sun")

    sim.gravity = "none" # Total isolation of the force
 
    rng = np.random.default_rng(42)
 
    occulter_hashes = []
    list_k = []
    for k in range(n_occulters):
        h = f"occ{k}"
        ang_k = 2*np.pi*k/n_occulters  
        x = AU*np.cos(ang_k)
        y = AU*np.sin(ang_k)
        sim.add(m=Mearth, x=x, y=y, z=0, hash=h)
        occulter_hashes.append(h)
        list_k.append(k)
 
    for i in range(N):
        x = AU + rng.uniform(-5, 5)*Rearth
        y = rng.uniform(-50, 50)*Rearth
        z = rng.uniform(-50, 50)*Rearth
        sim.add(m=0.0, x=x, y=y, z=z, vx=0, vy=0, vz=0)
 
    rebx = reboundx.Extras(sim)
    rad = rebx.load_force("radiation_forces")
    rebx.add_force(rad)
    rad.params["c"] = c_light
 
    sim.particles["sun"].params["radiation_source"] = 1
    sim.particles["sun"].params["R_eq"] = Rsun
 
    for h in occulter_hashes:
        sim.particles[h].params["shadow_creator"] = 1
        sim.particles[h].params["R_eq"] = Rearth*(0.5 + 0.2*list_k[occulter_hashes.index(h)])  
 
    for p in sim.particles[1+n_occulters:]:
        p.params["beta"] = 0.1
 
    if shadow_model != 0:
        rad.params["shadow_model"] = shadow_model

    sim._rebx_anchors = [rebx, rad]
    
    sim.dt = 1.0
    sim.integrator = "leapfrog"  
    return sim


def measure_performance(N, shadow_model, n_occulters=1, n_reps=7, steps_per_rep=20):
    sim = build_bench_sim(N, shadow_model, n_occulters)
    
    for _ in range(10):
        sim.step()
        
    times = []
    for _ in range(n_reps):
        t0 = time.perf_counter()
        for _ in range(steps_per_rep):
            sim.step()
        t1 = time.perf_counter()
        times.append((t1 - t0) / steps_per_rep)
        
    median_time = np.median(times)
    mad_time = np.median(np.abs(times - median_time))
    return median_time * 1000, mad_time * 1000  # ms



print("="*78)
print("BENCHMARK: average time per sim.step() [ms], ")
print("="*78)
print(f"{'N':>7} | {'Modelo':>12} | {'Mediana (ms)':>15} | {'MAD (ms)':>12} | {'Ratio vs None':>15}")
print("-" * 85)

for N in [10, 100, 1000, 5000, 20000]:

    models = [(0, "No Shadow"), (1, "Cylindrical"), (2, "Conical")]
    
    random.shuffle(models)

    results = {}
    for flag, name in models:
        med, mad = measure_performance(N, flag, n_occulters=1, n_reps=9, steps_per_rep=30)
        results[flag] = (name, med, mad)
        
    base_med = results[0][1]
    for flag in [0, 1, 2]:
        name, med, mad = results[flag]
        ratio = med / base_med
        print(f"{N:7d} | {name:>12} | {med:15.4f} | {mad:12.4f} | {ratio:15.2f}")
    print("-" * 85)

print("="*78)
print("BENCHMARK: escalating with nº occulting bodies (N=5000 fixed particles)")
print("="*78)
print(f"{'n_occulters':>12} | {'cylindrical (ms)':>16} | {'conical (ms)':>12}")
 
for n_occ in [1, 2, 4, 8]:
    med_cyl, _ = measure_performance(5000, 1, n_occulters=n_occ, n_reps=5, steps_per_rep=20)
    med_con, _ = measure_performance(5000, 2, n_occulters=n_occ, n_reps=5, steps_per_rep=20)
    print(f"{n_occ:12d} | {med_cyl:18.4f} | {med_con:15.4f}")