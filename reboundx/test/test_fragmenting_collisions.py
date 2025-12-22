import rebound
import reboundx
import unittest

class TestCollisionResolve(unittest.TestCase):
    def test_merging_collisions(self):
        # Merging
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
        com_f = sim.com()
        self.assertEqual(sim.N, 1)
        self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
        self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
        self.assertEqual(com_f.vy,0.0)
        self.assertEqual(com_f.vz,0.0)


    def test_elastic_bounce_collisions(self):
        # Elastic bounce
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
        self.assertEqual(sim.N, 2)
        com_f = sim.com()
        self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
        self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vy - com_f.vy)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vz - com_f.vz)/com_i.vx), 1e-10)


    def test_accreting_collisions(self):
        # Partial accretion
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
        print(sim.N)
        com_f = sim.com()
        self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
        self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vy - com_f.vy)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vz - com_f.vz)/com_i.vx), 1e-10)


    def test_erosive_collisions(self):
        # Non-grazing erosion
        sim = rebound.Simulation()
        rebx = reboundx.Extras(sim)
        collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
        rebx.add_collision_resolve(collision_resolve)
        collision_resolve.params["fc_min_frag_mass"] = 0.001
        sim.add(m=0.15, r=1.0, x=0)
        sim.add(m=0.10, r=1.0, x=10.0, vx=-30.0, vy=0.001, vz=0.001)
        sim.dt = 1
        sim.integrator = "mercurius"
        sim.collision = "direct"
        com_i = sim.com()
        self.assertEqual(sim.N, 2)
        sim.step()
        com_f = sim.com()
        self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
        self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vy - com_f.vy)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vz - com_f.vz)/com_i.vx), 1e-10)

        
    def test_grazing_collisions(self):
        # Grazing erosion
        sim = rebound.Simulation()
        rebx = reboundx.Extras(sim)
        collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
        rebx.add_collision_resolve(collision_resolve)
        collision_resolve.params["fc_min_frag_mass"] = 0.001
        sim.add(m=0.25, r=100.0, x=0)
        sim.add(m=0.20, r=100.0, x=200.0, y=150.0, vx=-100.0, vy=0.001, vz=0.001)
        sim.dt = 1
        sim.integrator = "mercurius"
        sim.collision = "direct"
        com_i = sim.com()
        self.assertEqual(sim.N, 2)
        sim.step()
        com_f = sim.com()
        self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
        self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vy - com_f.vy)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vz - com_f.vz)/com_i.vx), 1e-10)


    def test_erosive_collisions_integrators(self):
        # Non-grazing erosion, different integrators
        sim = rebound.Simulation()
        rebx = reboundx.Extras(sim)
        collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
        rebx.add_collision_resolve(collision_resolve)
        collision_resolve.params["fc_min_frag_mass"] = 0.001
        sim.add(m=0.15, r=1.0, x=0)
        sim.add(m=0.10, r=1.0, x=10.0, vx=-30.0, vy=0.001, vz=0.001)
        sim.dt = 1
        sim.collision = "direct"
        self.assertEqual(sim.N, 2)
        com_i = sim.com()

        integrators = ["mercurius", "ias15", "whfast", "BS"]
        #integrators = ["mercurius", "ias15", "whfast", "BS", "trace"]
        for integrator in integrators:
            sim.integrator = integrator
            sim.step()
            com_f = sim.com()
            self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
            self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
            self.assertLess(abs((com_i.vy - com_f.vy)/com_i.vx), 1e-10)
            self.assertLess(abs((com_i.vz - com_f.vz)/com_i.vx), 1e-10)
            print(f"integrator {integrator} passed.\n")
            

    def test_collision_with_the_sun(self):
        # Collision with the sun
        sim = rebound.Simulation()
        rebx = reboundx.Extras(sim)
        collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
        rebx.add_collision_resolve(collision_resolve)
        collision_resolve.params["fc_min_frag_mass"] = 0.001
        sim.add(m=1.0, r=0.1)
        sim.add(m=1e-7, r=1e-5, x=1.0, vx=-2.0, vy=0.001, vz=0.001)
        sim.dt = 1
        sim.integrator = "mercurius"
        sim.collision = "direct"
        com_i = sim.com()
        self.assertEqual(sim.N, 2)
        sim.step()
        com_f = sim.com()
        self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
        self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vy - com_f.vy)/com_i.vx), 1e-10)
        self.assertLess(abs((com_i.vz - com_f.vz)/com_i.vx), 1e-10)


    """
        def test_zero_mass(self):
            # One object with zero mass
            sim = rebound.Simulation()
            rebx = reboundx.Extras(sim)
            collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
            rebx.add_collision_resolve(collision_resolve)
            collision_resolve.params["fc_min_frag_mass"] = 0.001
            sim.add(m=1.0, r=0.1)
            sim.add(m=0.0, r=0.1, x=1.0, vx=-2.0, vy=0.001, vz=0.001)
            sim.dt = 1
            sim.integrator = "mercurius"
            sim.collision = "direct"
            com_i = sim.com()
            self.assertEqual(sim.N, 2)
            sim.step()
            print("N after zero = ", sim.N)
            com_f = sim.com()
            self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
            self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
            self.assertLess(abs((com_i.vy - com_f.vy)/com_i.vx), 1e-10)
            self.assertLess(abs((com_i.vz - com_f.vz)/com_i.vx), 1e-10)

        def test_no_min_frag_mass(self):
            # No min frag mass seg
            sim = rebound.Simulation()
            rebx = reboundx.Extras(sim)
            collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
            rebx.add_collision_resolve(collision_resolve)
            sim.add(m=0.15, r=1.0, x=0)
            sim.add(m=0.10, r=1.0, x=10.0, vx=-30.0, vy=0.001, vz=0.001)
            sim.dt = 1
            sim.integrator = "mercurius"
            sim.collision = "direct"
            com_i = sim.com()
            self.assertEqual(sim.N, 2)
            sim.step()
            print("N after no min frag mass = ", sim.N)
            com_f = sim.com()
            self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
            self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
            self.assertLess(abs((com_i.vy - com_f.vy)/com_i.vx), 1e-10)
            self.assertLess(abs((com_i.vz - com_f.vz)/com_i.vx), 1e-10)

        def test_zero_min_frag_mass(self):
            # Min frag mass = 0
            sim = rebound.Simulation()
            rebx = reboundx.Extras(sim)
            collision_resolve = rebx.load_collision_resolve("fragmenting_collisions")
            rebx.add_collision_resolve(collision_resolve)
            collision_resolve.params["fc_min_frag_mass"] = 0.0
            sim.add(m=0.15, r=1.0, x=0)
            sim.add(m=0.10, r=1.0, x=10.0, vx=-30.0, vy=0.001, vz=0.001)
            sim.dt = 1
            sim.integrator = "mercurius"
            sim.collision = "direct"
            com_i = sim.com()
            self.assertEqual(sim.N, 2)
            sim.step()
            print("N after no min frag mass = ", sim.N)
            com_f = sim.com()
            self.assertLess(abs((com_i.m - com_f.m)/com_i.m), 1e-10)
            self.assertLess(abs((com_i.vx - com_f.vx)/com_i.vx), 1e-10)
            self.assertLess(abs((com_i.vy - com_f.vy)/com_i.vx), 1e-10)
            self.assertLess(abs((com_i.vz - com_f.vz)/com_i.vx), 1e-10)
    """



if __name__ == '__main__':
    unittest.main()
