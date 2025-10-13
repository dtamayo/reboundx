import rebound
import reboundx
import unittest

class TestCollisionResolve(unittest.TestCase):

    def test_merging_collisions(self):
        sim = rebound.Simulation()
        rebx = reboundx.Extras(sim)
        collision_resolve = rebx.load_collision_resolve("merging_collisions")
        rebx.add_collision_resolve(collision_resolve)
        sim.add(r=1)
        sim.add(r=1,x=2.5,vx=-1)
        sim.dt = 1
        sim.integrator = "leapfrog"
        sim.collision = "direct"
        self.assertEqual(sim.N, 2)
        sim.step()
        self.assertEqual(sim.N, 1)


if __name__ == '__main__':
    unittest.main()
