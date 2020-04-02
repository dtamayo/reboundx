import rebound
import reboundx
import unittest
import os
import numpy as np

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
binary = os.path.join(THIS_DIR, 'binaries/twoplanets.bin')

class TestConservation(unittest.TestCase):
    def test_adiabatic_mass_loss(self):
        sim = rebound.Simulation(binary)
        sim.integrator = "ias15"
        rebx = reboundx.Extras(sim)
        times = [0, 2000., 4000., 6000., 8000., 10000.]
        values = [1., 0.8, 0.6, 0.4, 0.3, 0.2]
        starmass = reboundx.Interpolator(rebx, times, values, "spline")

        Nout = 1000
        ts = np.linspace(0., 1.e4, Nout)
        ps  = sim.particles
        a10 = ps[1].a
        for i, time in enumerate(ts):
            sim.integrate(time)
            ps[0].m = starmass.interpolate(rebx, t=sim.t)
            sim.move_to_com() # lost mass had momentum, so need to move back to COM frame
         
        self.assertLess(abs((ps[0].m-values[-1])/values[-1]), 1.e-6)
        self.assertLess(abs((ps[1].a-5*a10)/a10), 1.e-2)
   
    def test_klo(self):
        sim = rebound.Simulation(binary)
        sim.integrator = "ias15"
        rebx = reboundx.Extras(sim)
        times = [0, 2000., 4000., 6000., 8000., 10000.]
        values = [1., 0.8, 0.6, 0.4, 0.3, 0.2]
        starmass = reboundx.Interpolator(rebx, times, values, "spline")
        
        m0 = sim.particles[0].m
        mint0 = starmass.interpolate(rebx, t=0)

        self.assertLess(abs((m0-mint0)/m0), 1.e-6)

        mint1 = starmass.interpolate(rebx, t=5000) # advance klo
        mint1 = starmass.interpolate(rebx, t=0) # advance klo
        
        self.assertLess(abs((m0-mint1)/m0), 1.e-6)

    def test_backwards_time(self):
        sim = rebound.Simulation(binary)
        sim.integrator = "ias15"
        rebx = reboundx.Extras(sim)
        times = [0, 2000., 4000., 6000., 8000., 10000.]
        values = [1., 0.8, 0.6, 0.4, 0.3, 0.2]
        starmass = reboundx.Interpolator(rebx, times, values, "spline")

        Nout = 1000
        ts = np.linspace(0., 1.e4, Nout)
        ps  = sim.particles
        a10 = ps[1].a
        m0 = ps[0].m

        for i, time in enumerate(ts):
            sim.integrate(time)
            ps[0].m = starmass.interpolate(rebx, t=sim.t)
            sim.move_to_com() # lost mass had momentum, so need to move back to COM frame
        sim.dt = -0.001
        for i, time in enumerate(times[::-1]):
            sim.integrate(time)
            ps[0].m = starmass.interpolate(rebx, t=sim.t)
            sim.move_to_com() # lost mass had momentum, so need to move back to COM frame
        self.assertLess(abs((ps[0].m-m0)/m0), 1.e-2)
        self.assertLess(abs((ps[1].a-a10)/a10), 1.e-2)
    
if __name__ == '__main__':
    unittest.main()
