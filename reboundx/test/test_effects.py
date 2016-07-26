import rebound
import reboundx
import unittest
import os

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
        sources = [each for each in os.listdir('../../src/') if each.endswith('.c')] 
        for source in sources:
            if source == "core.c" or source == "rebxtools.c":
                continue
            effect = self.rebx.add(source[:-2])

if __name__ == '__main__':
    unittest.main()
