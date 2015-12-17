import rebound
import reboundx
from reboundx import data
import unittest
import math

# when we add pericenter node precession testing, keep in mind the abs. error should be ~mass ratio/e (or i) due to 1/ei in Lagrange eq.
class Test_modify_orbits_forces(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        data.add_earths(self.sim, ei=1.e-3) # ei=val of e&i. Forces wrong at order ei^2, so choose low enough that 
                                            # it doesn't matter at the level of the planet/star mass ratio
        self.sim.move_to_com()
        earth = self.sim.particles[1]
        earthP0 = earth.P
        earthPf = earthP0/math.e**1.5
        self.sim.dt = 0.05123*earthPf
        self.tmax = 50.*earthP0             # need to choose short tmax compared to t_sec (see research notebook)
        self.rebx = reboundx.Extras(self.sim)
        self.abstolfac = 2.*math.pi # factor to multiply order of magnitude abs err estimate
    
    def tearDown(self):
        self.sim = None
    
    def test_migration_ias15(self):
        """
        Set tau_a to -tmax and check that a at the end is a0/e to within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "ias15"
        self.rebx.add_modify_orbits_forces()
        earth = self.sim.particles[1]
        earth.tau_a = -self.tmax
        a0 = earth.a
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.a - a0/math.e), earth.m/sun.m*self.abstolfac)
        # osculating elements wrong at level of planet/star mass ratio (for separation/a ~1), so this sets absolute tolerance.
    
    def test_migration_whfast(self):
        """
        Set tau_a to -tmax and check that a at the end is a0/e to within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "whfast"
        self.rebx.add_modify_orbits_forces()
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
        self.rebx.add_modify_orbits_forces()
        earth = self.sim.particles[1]
        earth.tau_e = -self.tmax
        e0 = earth.e
        eExp = e0/math.e
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.e - e0/math.e), earth.m/sun.m*self.abstolfac)
    
    def test_eDamping_whfast(self):
        """
        Set tau_e to -tmax and check that e at the end is e0/e to within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "whfast"
        self.rebx.add_modify_orbits_forces()
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
        self.rebx.add_modify_orbits_forces()
        earth = self.sim.particles[1]
        earth.tau_inc = -self.tmax
        i0 = earth.inc
        eExp = i0/math.e
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.inc - i0/math.e), earth.m/sun.m*self.abstolfac)
    
    def test_incDamping_whfast(self):
        """
        Set tau_inc to -tmax and check that i at the end is i0/e to within ~ mass ratio (3e-6)
        """
        self.sim.integrator = "whfast"
        self.rebx.add_modify_orbits_forces()
        earth = self.sim.particles[1]
        earth.tau_inc = -self.tmax
        i0 = earth.inc
        self.sim.integrate(self.tmax)
        sun = self.sim.particles[0]
        self.assertLess(abs(earth.inc - i0/math.e), earth.m/sun.m*self.abstolfac)
    
if __name__ == '__main__':
    unittest.main()
