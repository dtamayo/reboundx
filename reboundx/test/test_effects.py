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
        rf = self.rebx.add_radiation_forces()
        mod = self.rebx.add_modify_orbits_direct()
        mof = self.rebx.add_modify_orbits_forces()
        gr = self.rebx.add_gr()
        grf = self.rebx.add_gr_full()
        grp = self.rebx.add_gr_potential()

        rf.c = 3.
        self.assertAlmostEqual(rf.c, 3.)
        rf.source_index = 2
        self.assertEqual(rf.source_index, 2)

        gr.c = 3.
        self.assertAlmostEqual(gr.c, 3.)
        gr.source_index = 2
        self.assertEqual(gr.source_index, 2)
        
        grp.c = 3.
        self.assertAlmostEqual(grp.c, 3.)
        grp.source_index = 2
        self.assertEqual(grp.source_index, 2)
        
        grf.c = 3.
        self.assertAlmostEqual(grf.c, 3.)

        mod.p = 0.2
        self.assertAlmostEqual(mod.p, 0.2)
        mod.coordinates = 3
        self.assertEqual(mod.coordinates, 3)

        mof.coordinates = 3
        self.assertEqual(mof.coordinates, 3)

if __name__ == '__main__':
    unittest.main()
