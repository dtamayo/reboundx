import rebound
import reboundx
import unittest

class TestCollisionResolve(unittest.TestCase):
    def test_merging_collisions(self):
        print("test merge started")
        sim = rebound.Simulation()
        rebx = reboundx.Extras(sim)
        collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
        rebx.add_collision_resolve(collision_resolve)
        collision_resolve.params["fc_min_frag_mass"] = 1.0
        sim.add(m=1, r=1)
        sim.add(m=1, r=1,x=2.5,vx=-1)
        sim.dt = 1
        sim.integrator = "mercurius"
        sim.collision = "direct"
        com_i = sim.com()
        self.assertEqual(sim.N, 2)
        sim.step()
        print('test merge integrated')
        com_f = sim.com()
        self.assertEqual(sim.N, 1)
        self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
        self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
        self.assertEqual(com_f.vy,0.0)
        self.assertEqual(com_f.vz,0.0)


    def test_elastic_bounce_collisions(self):
        print('test elastic bounce started')
        sim = rebound.Simulation()
        rebx = reboundx.Extras(sim)
        collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
        rebx.add_collision_resolve(collision_resolve)
        collision_resolve.params["fc_min_frag_mass"] = 0.05
        sim.add(m=0.101, r=0.5)
        sim.add(m=0.10, r=0.5, x=1.0, y=0.9, vx=-1.0, vy=0.001, vz=0.001)   

        sim.dt = 1
        sim.integrator = "mercurius"
        sim.collision = "direct"
        com_i = sim.com()
        self.assertEqual(sim.N, 2)
        sim.step()
        print('test eb integrated')
        self.assertEqual(sim.N, 2)
        com_f = sim.com()
        self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
        self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vy - com_f.vy)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vz - com_f.vz)/com_i.vx), 1e-10)


    def test_fragmenting_collisions(self):
        print("test frag started")
        sim = rebound.Simulation()
        rebx = reboundx.Extras(sim)
        collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
        rebx.add_collision_resolve(collision_resolve)
        collision_resolve.params["fc_min_frag_mass"] = 0.001
        sim.add(m=0.15, r=1.0, x=0)
        sim.add(m=0.10, r=1.0, x=2.0, vx=-3.0, vy=0.001, vz=0.001)
        sim.dt = 1
        sim.integrator = "mercurius"
        sim.collision = "direct"
        com_i = sim.com()
        self.assertEqual(sim.N, 2)
        sim.step()
        print('test frag integrated')
        print(sim.N)
        #self.assertGreater(sim.N, 2)
        com_f = sim.com()
        self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
        self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vy - com_f.vy)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vz - com_f.vz)/com_i.vx), 1e-10)




if __name__ == '__main__':
    unittest.main()
