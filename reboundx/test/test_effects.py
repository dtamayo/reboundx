import rebound
import reboundx
import unittest

class Test_effects(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        rebx = reboundx.Extras(self.sim)

    def tearDown(self):
        self.sim = None
    
    def test_effects(self):
        """
        Assign all the particle to two particles, and make sure we get the same values back
        """
        rebx.add_radiation_forces()
        rebx.add_modify_orbits_direct()
        rebx.add_modify_orbits_forces()
        rebx.add_gr()
        rebx.add_gr_full()
        rebx.add_gr_potential()

        rebx.radiation_forces.c = 3.
        self.assertAlmostEqual(rebx.radiation_forces.c, 3.)
        rebx.radiation_forces.source_index = 2
        self.assertEqual(rebx.radiation_forces.source_index, 2)

        rebx.gr.c = 3.
        self.assertAlmostEqual(rebx.gr.c, 3.)
        rebx.gr.source_index = 2
        self.assertEqual(rebx.gr.source_index, 2)
        
        rebx.gr_potential.c = 3.
        self.assertAlmostEqual(rebx.gr_potential.c, 3.)
        rebx.gr_potential.source_index = 2
        self.assertEqual(rebx.gr_potential.source_index, 2)
        
        rebx.gr_full.c = 3.
        self.assertAlmostEqual(rebx.gr_full.c, 3.)

        rebx.modify_orbits_direct.p = 0.2
        self.assertAlmostEqual(rebx.modify_orbits_direct.p, 0.2)
        rebx.modify_orbits_direct.coordinates = 3
        self.assertEqual(rebx.modify_orbits_direct.coordinates, 3)

        rebx.modify_orbits_forces.p = 0.2
        self.assertAlmostEqual(rebx.modify_orbits_forces.p, 0.2)
        rebx.modify_orbits_forces.coordinates = 3
        self.assertEqual(rebx.modify_orbits_forces.coordinates, 3)

if __name__ == '__main__':
    unittest.main()
