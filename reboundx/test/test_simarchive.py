import rebound
import reboundx
import unittest

class TestReproducibility(unittest.TestCase):
    def test_reproducibility(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        rebx = reboundx.Extras(sim)
        gr = rebx.add("gr")
        gr.params['c'] = 1.e3
        sim.automateSimulationArchive("test.bin",interval=1e3,deletefile=True)
        sim.integrate(1.e4)
        rebx.save("rebx.bin")	
        sa = rebound.SimulationArchive("test.bin", rebxfilename="rebx.bin")

        for i in [2,5,7]:
            sim = sa[i]
            x = sim.particles[1].x
            sim = sa[i-1]
            sim.integrate(sa[i].t, exact_finish_time=0)
            self.assertEqual(x, sim.particles[1].x)

if __name__ == '__main__':
    unittest.main()
