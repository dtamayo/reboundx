import rebound
import reboundx
import unittest

class TestCollisionResolve(unittest.TestCase):

    def test_merging_collisions(self):
        sim = rebound.Simulation()
        rebx = reboundx.Extras(sim)
        collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
        rebx.add_collision_resolve(collision_resolve)
        sim.add(r=1)
        sim.add(r=1,x=2.5,vx=-1)
        sim.dt = 1
        sim.integrator = "leapfrog"
        sim.collision = "direct"
        com_i = sim.com()
        self.assertEqual(sim.N, 2)
        sim.step()
        com_f = sim.com()
        self.assertEqual(sim.N, 1)
        self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
        self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vy - com_f.vy)/com_i.vy), 1e-10)
        self.assertLess(abs((com_i.vz - com_f.vz)/com_i.vz), 1e-10)




if __name__ == '__main__':
    unittest.main()
