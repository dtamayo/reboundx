import rebound
import reboundx
from reboundx import data
import unittest
import math

class Test_modify_orbits_direct(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        data.add_earths(self.sim, ei=0.1) # Orbits direct exact, need to set e&i high enough that omega/Omega don't fly around
        self.sim.move_to_com()
        self.integrators = ["ias15", "whfast"]
        earth = self.sim.particles[1]
        earthP0 = earth.P
        earthPf = earthP0/math.e**1.5
        self.sim.dt = 0.05123*earthPf
        self.tmax = 50.*earthP0
        self.rebx = reboundx.Extras(self.sim)
        self.abstolfac = 2*math.pi # factor to multiply order of magnitude abs err estimate
    
    def tearDown(self):
        self.sim = None
    
    def test_migration_ias15(self):
        """
        Set tau_a to -tmax and check that a at the end is a0/e to within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "ias15"
        self.rebx.add_modify_orbits_direct()
        earth = self.sim.particles[1]
        earth.tau_a = -self.tmax
        a0 = earth.a
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.a - a0/math.e), earth.m/sun.m*self.abstolfac)
    
    def test_migration_whfast(self):
        """
        Set tau_a to -tmax and check that a at the end is a0/e to within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "whfast"
        self.rebx.add_modify_orbits_direct()
        earth = self.sim.particles[1]
        earth.tau_a = -self.tmax
        a0 = earth.a
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.a - a0/math.e), earth.m/sun.m*self.abstolfac)

    def test_eDamping_ias15(self):
        """
        Set tau_e to -tmax and check that e at the end is e0/e to within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "ias15"
        self.rebx.add_modify_orbits_direct()
        earth = self.sim.particles[1]
        earth.tau_e = -self.tmax
        e0 = earth.e
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.e - e0/math.e), earth.m/sun.m*self.abstolfac)
    
    def test_eDamping_whfast(self):
        """
        Set tau_e to -tmax and check that e at the end is e0/e to within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "whfast"
        self.rebx.add_modify_orbits_direct()
        earth = self.sim.particles[1]
        earth.tau_e = -self.tmax
        e0 = earth.e
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.e - e0/math.e), earth.m/sun.m*self.abstolfac)

    def test_incDamping_ias15(self):
        """
        Set tau_inc to -tmax and check that i at the end is i0/e to within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "ias15"
        self.rebx.add_modify_orbits_direct()
        earth = self.sim.particles[1]
        earth.tau_inc = -self.tmax
        i0 = earth.inc
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.inc - i0/math.e), earth.m/sun.m*self.abstolfac)
    
    def test_incDamping_whfast(self):
        """
        Set tau_inc to -tmax and check that i at the end is i0/e to within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "whfast"
        self.rebx.add_modify_orbits_direct()
        earth = self.sim.particles[1]
        earth.tau_inc = -self.tmax
        i0 = earth.inc
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.inc - i0/math.e), earth.m/sun.m*self.abstolfac)

    def test_omegaPrec_ias15(self):
        """
        Set tau_omega to tmax and check that it returns to omega=omega0  at the end within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "ias15"
        self.rebx.add_modify_orbits_direct()
        earth = self.sim.particles[1]
        earth.tau_omega = self.tmax
        omega0 = earth.omega
        y0 = earth.y
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.omega - omega0), earth.m/sun.m*self.abstolfac/earth.e)
    
    def test_omegaPrec_whfast(self):
        """
        Set tau_omega to tmax and check that it returns to omega=omega0  at the end within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "whfast"
        self.rebx.add_modify_orbits_direct()
        earth = self.sim.particles[1]
        earth.tau_omega = self.tmax
        omega0 = earth.omega
        a0 = earth.a
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.omega - omega0), earth.m/sun.m*self.abstolfac/earth.e)
    
    def test_OmegaPrec_ias15(self):
        """
        Set tau_Omega to tmax and check that it returns to Omega=Omega0  at the end within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "ias15"
        self.rebx.add_modify_orbits_direct()
        earth = self.sim.particles[1]
        earth.tau_Omega = self.tmax
        Omega0 = earth.Omega
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.Omega - Omega0), earth.m/sun.m*self.abstolfac/earth.inc)
    
    def test_OmegaPrec_whfast(self):
        """
        Set tau_Omega to tmax and check that it returns to Omega=Omega0  at the end within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "whfast"
        self.rebx.add_modify_orbits_direct()
        earth = self.sim.particles[1]
        earth.tau_Omega = self.tmax
        Omega0 = earth.Omega
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.Omega - Omega0), earth.m/sun.m*self.abstolfac/earth.inc)

if __name__ == '__main__':
    unittest.main()
