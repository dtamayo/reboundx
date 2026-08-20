import rebound
import reboundx
import numpy as np
import matplotlib.pyplot as plt

m_sun = 1.989e30
m_earth = 5.972e24
m_mercuri = 3.301e23
Rsun = 6.957e8
Rearth = 6.371e6
Rmercuri = 2.439e6
c_light = 3e8  
r0 = 7000e3  
AU = 1.496e11

d_mercuri = 5.79e10 


def build_simulation(shadow_flag, eps):  

    sim = rebound.Simulation()
    sim.units = ('m', 's', 'kg')
    rebx = reboundx.Extras(sim)

    rf = rebx.load_force("radiation_forces")
    rf.params["c"] = c_light  
    rf.params["shadow_model"] = shadow_flag  # 0 = No shadow, 1 = Cylindrical, 2 = Conical

    # 0. Light Source (Sun)
    sim.add(m=m_sun, x=0.0, y=0.0, z=0.0, hash="Sun")
    Sun = sim.particles["Sun"]
    Sun.params["radiation_source"] = 1
    Sun.params["R_eq"] = Rsun

    # 1. Occulting body (Earth)
    v_earth = np.sqrt(sim.G * m_sun / AU)

    sim.add(m=m_earth, x=AU, y=0.0, z=0.0, vx=0.0, vy=v_earth, vz=0.0, hash="Earth")
    Earth = sim.particles["Earth"]
    Earth.params["shadow_creator"] = 1
    #Earth.params["R_eq"] = Rearth

    # 1. (Mercury)
    v_mercury = np.sqrt(sim.G * m_sun / d_mercuri)
    sim.add(m=m_mercuri, x=d_mercuri, y=0.0, z=0.0, vx=0.0, vy=v_mercury, vz=0.0, hash="Mercury")
    Mercury = sim.particles["Mercury"]
    Mercury.params["shadow_creator"] = 1
    Mercury.params["R_eq"] = Rmercuri

    # 2. Satellite
    mu = sim.G * m_earth
    v0_sat = np.sqrt(mu / r0)
    sim.add(m=1000.0, x=r0 + AU, y=0.0, z=0.0, vx=Earth.vx, vy=Earth.vy + v0_sat, vz=Earth.vz, hash= "Sat")
    Sat = sim.particles["Sat"]
    Sat.params["beta"] = 1e-4

    rebx.add_force(rf)

    sim.integrator = "ias15"
    sim.ri_ias15.epsilon = eps
    sim.rebx_anchor = [rebx, rf]
    return sim

def get_shadow_factor(sim, shadow_model):
    """
    Calculates if the satellite is in geometric shadow evaluating ALL occulting bodies.
    If the shadow_model is conical (2), returns True when it touches gloom.
    """
    if shadow_model == 0:
        return 1.0
            
    sun = sim.particles["Sun"]
    sat = sim.particles["Sat"]
    R_source = sun.params.get("R_eq", 0.0)
    
    rx = sat.x - sun.x
    ry = sat.y - sun.y
    rz = sat.z - sun.z
    r_mag = np.sqrt(rx**2 + ry**2 + rz**2)
    
    shadow_factor = 1.0
    
    for p in sim.particles:
        if p.index == sun.index or p.index == sat.index:
            continue
            
        if p.params.get("shadow_creator") != 1:
            continue
            
        try:
            R_occ = p.params["R_eq"]
        except AttributeError:
            continue
            
        dx_occ = p.x - sun.x
        dy_occ = p.y - sun.y
        dz_occ = p.z - sun.z
        dr_occ = np.sqrt(dx_occ**2 + dy_occ**2 + dz_occ**2)
        
        if dr_occ <= 0:
            continue
        
        if shadow_model == 1:
            ux = dx_occ / dr_occ
            uy = dy_occ / dr_occ
            uz = dz_occ / dr_occ
            s = rx*ux + ry*uy + rz*uz
            
            if s <= dr_occ:
                continue 
            
            perp = np.sqrt((rx - s*ux)**2 + (ry - s*uy)**2 + (rz - s*uz)**2)
            if perp < R_occ:
                shadow_factor = 0.0
                break 
                
        elif shadow_model == 2:
            r_px = sat.x - p.x
            r_py = sat.y - p.y
            r_pz = sat.z - p.z
            r_pmod = np.sqrt(r_px**2 + r_py**2 + r_pz**2)
            
            dot = -(dx_occ*r_px + dy_occ*r_py + dz_occ*r_pz)
            s0 = -dot / dr_occ
            
            if s0 < 0.0:
                continue
            
            l_arg = max(0.0, r_pmod**2 - s0**2)
            l = np.sqrt(l_arg)
            
            app_r_source = np.arcsin(min(R_source / r_mag, 1.0))
            app_r_occ = np.arcsin(min(R_occ / r_pmod, 1.0))
            d_ang = np.arcsin(min(l / r_pmod, 1.0))
            
            f_local = 1.0
            if d_ang >= app_r_source + app_r_occ:
                f_local = 1.0
            elif d_ang <= abs(app_r_source - app_r_occ):
                if app_r_source <= app_r_occ: 
                    f_local = 0.0
                else: 
                    f_local = 1.0 - (app_r_occ**2) / (app_r_source**2)
            else:
                d1 = (app_r_source**2 - app_r_occ**2 + d_ang**2) / (2.0 * d_ang)
                d2 = d_ang - d1
                arg1 = max(-1.0, min(1.0, d1 / app_r_source))
                arg2 = max(-1.0, min(1.0, d2 / app_r_occ))
                
                area = (app_r_source**2 * np.arccos(arg1) - d1 * np.sqrt(abs(app_r_source**2 - d1**2)) +
                        app_r_occ**2 * np.arccos(arg2) - d2 * np.sqrt(abs(app_r_occ**2 - d2**2)))
                f_local = 1.0 - area / (np.pi * app_r_source**2)
            
            if f_local < shadow_factor:
                shadow_factor = f_local

    return shadow_factor


n_orbits = 50
G = 6.674e-11
T = 2*np.pi*np.sqrt(r0**3 / (G*m_earth))
t_end = n_orbits * T
epsilons = [1e-7, 1e-9, 1e-11]

results = {}

print("="*80)
print("Convergence test and shadowing dynamics (IAS15)")
print("="*80)

for shadow_flag, name in [(0, "no shadow"), (1, "cylindrical"), (2, "conical")]:

    results[name] = {}
    print(f"\n[{name.upper()}]")
    print(f"{'Epsilon':>9} | {'N_steps':>7} | {'Delta a (m)':>18} | {'Eclipse Dur (s)':>15} | {'<f> sombra':>10}")
    print("-" * 80)

    for eps in epsilons: 

        sim = build_simulation(shadow_flag, eps)
        a_inicial = sim.particles["Sat"].orbit(primary=sim.particles["Earth"]).a

        ts, dts, a_vals = [], [], []
        n_steps = 0
        eclipse_time = 0.0
        sum_f_shadow = 0.0
        shadow_steps = 0

        while sim.t < t_end:
            sim.step()
            n_steps += 1

            ts.append(sim.t)
            dt = sim.dt
            dts.append(dt)
            orb = sim.particles["Sat"].orbit(primary=sim.particles["Earth"])
            a_vals.append(orb.a)
            f = get_shadow_factor(sim, shadow_flag)
            if f < 1.0:
                eclipse_time += dt
                sum_f_shadow += f * dt  
    
        ts = np.array(ts); dts = np.array(dts); a_vals = np.array(a_vals)
    

        avg_f = (sum_f_shadow / eclipse_time) if eclipse_time > 0 else 1.0

        a_final = sim.particles["Sat"].orbit(primary=sim.particles["Earth"]).a
        delta_a = a_final - a_inicial
    
        results[name][eps] = dict(t=ts, dt=dts, a=a_vals, delta_a= delta_a)
        print(f"{eps:9.1e} | {n_steps:7d} | {delta_a:18.10e} | {eclipse_time:15.2f}| {avg_f:10.5f}")


print("\n" + "="*80)
print("Convergence analysis (Epsilon = 1e-11)")
print("="*80)

da_none = results["no shadow"][1e-11]["delta_a"]
da_cyl  = results["cylindrical"][1e-11]["delta_a"]
da_con  = results["conical"][1e-11]["delta_a"]

print(f"Delta a (No Shadow)   : {da_none:18.10e} m")
print(f"Delta a (Cylindrical) : {da_cyl:18.10e} m")
print(f"Delta a (Conical)     : {da_con:18.10e} m")
print("-" * 80)
print(f"Cylindrical effect (None - Cyl) : {da_none - da_cyl:18.10e} m")
print(f"Conical effect (None - Con)     : {da_none - da_con:18.10e} m")
print(f"Difference between Cylindrical vs Conical: {da_cyl - da_con:18.10e} m")

print("\n" + "="*70)
print("COMPARISON OF RESULTS")
print("="*70)

a_none = results["no shadow"][1e-11]["a"][-1]
a_cyl  = results["cylindrical"][1e-11]["a"][-1]
a_con  = results["conical"][1e-11]["a"][-1]

print(f"delta_a (no shadow - cylindrical) : {(a_none - a_cyl):.6e} m")
print(f"delta_a (no shadow - conical)     : {(a_none - a_con):.6e} m")
print(f"delta_a (cylindrical - conical)   : {(a_cyl - a_con):.6e} m")
print(f"\ndt minimum in cylindrical shadow   : {results['cylindrical'][1e-11]['dt'].min():.4e} s")
print(f"dt minimum in conical shadow       : {results['conical'][1e-11]['dt'].min():.4e} s")



plt.figure()
plt.plot(results["no shadow"][1e-11]["t"], results["no shadow"][1e-11]["a"], label="No Shadow", color="black")
plt.plot(results["cylindrical"][1e-11]["t"], results["cylindrical"][1e-11]["a"], label="Cylindrical Shadow", color="blue")
plt.plot(results["conical"][1e-11]["t"], results["conical"][1e-11]["a"], label="Conical Shadow", color="red")
plt.xlabel("Time (s)")
plt.ylabel("Semi-major Axis (m)")
plt.title("Evolution of the semi-major axis of the satellite's orbit")
plt.legend()
plt.show()