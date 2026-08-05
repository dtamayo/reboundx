import rebound
import reboundx
import numpy as np

def test_radiation_shadows():
    """
    Test the radiation forces shadowing models (0: None, 1: Cylindrical, 2: Conical).
    Validates that particles in the umbra feel 0 acceleration, and particles 
    in full light feel the exact theoretical radiation acceleration.
    """
    sim = rebound.Simulation()
    sim.units = ('m', 's', 'kg')
    
    # Disable gravity to mathematically isolate the radiation force
    sim.gravity = "none"
    
    rebx = reboundx.Extras(sim)
    rf = rebx.load_force("radiation_forces")
    rf.params["c"] = 299792458.0
    rf.params["shadow_model"] = 2  # Activate Conical Model

    # Light Source (Sun)
    sim.add(m=1.989e30, x=0, y=0, z=0, hash="Sun")
    sim.particles["Sun"].params["radiation_source"] = 1
    sim.particles["Sun"].params["R_eq"] = 6.957e8

    # Occulting body (Earth) at 1 AU
    au = 1.496e11
    sim.add(m=5.972e24, x=au, y=0, z=0, hash="Earth")
    sim.particles["Earth"].params["shadow_creator"] = 1
    sim.particles["Earth"].params["R_eq"] = 6.371e6

    # Test Particle 1: Deep in the umbra (just 10,000 km behind Earth)
    sim.add(m=0, x=au + 10000e3, y=0, z=0, hash="Dust_Umbra")
    sim.particles["Dust_Umbra"].params["beta"] = 0.1

    # Test Particle 2: In full light (10,000 km in front of Earth)
    sim.add(m=0, x=au - 10000e3, y=0, z=0, hash="Dust_Light")
    sim.particles["Dust_Light"].params["beta"] = 0.1

    rebx.add_force(rf)
    
    sim.step()
    
    # 1. Verification of Full Shadow (Umbra)
    a_umbra = sim.particles["Dust_Umbra"].ax
    assert a_umbra == 0.0, f"Expected 0.0 acceleration in umbra, got {a_umbra}"
    
    # 2. Verification of Full Light
    # The acceleration must match the exact theoretical PR formula (a = beta * GM / r^2)
    dist_light = sim.particles["Dust_Light"].x
    expected_a_light = 0.1 * sim.G * sim.particles["Sun"].m / (dist_light**2)
    a_light = sim.particles["Dust_Light"].ax
    
    # Use np.isclose for floating point comparisons
    assert np.isclose(a_light, expected_a_light, rtol=1e-7), \
        f"Expected {expected_a_light} in full light, got {a_light}"