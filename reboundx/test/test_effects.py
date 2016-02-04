import rebound
import reboundx
import unittest

class Test_effects(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.rebx = reboundx.Extras(self.sim)

    def tearDown(self):
        self.sim = None
    
    def test_effects(self):
        """
        Assign all the particle to two particles, and make sure we get the same values back
        """
        self.rebx.add_radiation_forces()
        self.rebx.add_modify_orbits_direct()
        self.rebx.add_modify_orbits_forces()
        self.rebx.add_gr()
        self.rebx.add_gr_full()
        self.rebx.add_gr_potential()

        self.rebx.radiation_forces.c = 3.
        self.assertAlmostEqual(self.rebx.radiation_forces.c, 3.)
        self.rebx.radiation_forces.source_index = 2
        self.assertEqual(self.rebx.radiation_forces.source_index, 2)

        self.rebx.gr.c = 3.
        self.assertAlmostEqual(self.rebx.gr.c, 3.)
        self.rebx.gr.source_index = 2
        self.assertEqual(self.rebx.gr.source_index, 2)
        
        self.rebx.gr_potential.c = 3.
        self.assertAlmostEqual(self.rebx.gr_potential.c, 3.)
        self.rebx.gr_potential.source_index = 2
        self.assertEqual(self.rebx.gr_potential.source_index, 2)
        
        self.rebx.gr_full.c = 3.
        self.assertAlmostEqual(self.rebx.gr_full.c, 3.)

        self.rebx.modify_orbits_direct.p = 0.2
        self.assertAlmostEqual(self.rebx.modify_orbits_direct.p, 0.2)
        self.rebx.modify_orbits_direct.coordinates = 3
        self.assertEqual(self.rebx.modify_orbits_direct.coordinates, 3)

        self.rebx.modify_orbits_forces.p = 0.2
        self.assertAlmostEqual(self.rebx.modify_orbits_forces.p, 0.2)
        self.rebx.modify_orbits_forces.coordinates = 3
        self.assertEqual(self.rebx.modify_orbits_forces.coordinates, 3)

if __name__ == '__main__':
    unittest.main()
