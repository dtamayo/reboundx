import rebound
import reboundx
import unittest
import numpy as np

TMAX = 1.e3
TAU_INNER = -TMAX       # don't want too many e-foldings so e not too small
TAU_OUTER = -2. * TMAX
NOUT = 1000             # discrepancies between dt and dt_last_done largest at
                        # exact finish times, so worst case is stopping many times

def make_sim(integrator):
    sim = rebound.Simulation()
    sim.add(m=1.)
    sim.add(m=1e-6, a=1, e=0.1, inc=0.1)
    sim.add(m=1e-6, a=10, e=0.1, inc=0.1)
    sim.dt = sim.particles[1].P / 123.0
    sim.integrator = integrator
    sim.move_to_com()
    return sim

def add_damping(sim):
    """Attach modify_orbits_direct. Returns the rebx instance (kept alive by
    the caller, since Extras holds the only reference to the operator)."""
    rebx = reboundx.Extras(sim)
    mod = rebx.load_operator("modify_orbits_direct")
    rebx.add_operator(mod)
    for p, tau in ((sim.particles[1], TAU_INNER), (sim.particles[2], TAU_OUTER)):
        p.params["tau_e"] = tau
        p.params["tau_inc"] = tau
    return rebx

class TestModifyOrbitDirect(unittest.TestCase):

    def check_damped(self, integrator):
        """Elements must follow x = x0*exp(t/tau) to within 1%."""
        sim = make_sim(integrator)
        rebx = add_damping(sim)
        for t in np.linspace(0., TMAX, NOUT):
            sim.integrate(t)

        for i, tau in ((1, TAU_INNER), (2, TAU_OUTER)):
            for element, x0 in (("e", 0.1), ("inc", 0.1)):
                want = x0 * np.exp(TMAX / tau)
                got = getattr(sim.particles[i], element)
                with self.subTest(integrator=integrator, particle=i, element=element):
                    self.assertAlmostEqual(
                        got / want, 1.0, delta=1e-2,
                        msg="%s: particle %d %s = %.8e, expected %.8e "
                            "(relative error %.2e, tolerance %.1e)"
                            % (integrator, i, element, got, want,
                               abs(got / want - 1.), 1e-2))

    def check_unsupported(self, integrator):
        sim = make_sim(integrator)
        rebx = reboundx.Extras(sim)
        mod = rebx.load_operator("modify_orbits_direct")
        with self.assertRaises(RuntimeError) as cm:
            rebx.add_operator(mod)

    def check_noop(self, integrator):
        """Orbit must be untouched: no damping, and no motion either."""
        sim = make_sim(integrator)
        rebx = add_damping(sim)
        before = [(p.a, p.e, p.inc) for p in sim.particles[1:]]
        for t in np.linspace(0., TMAX, NOUT):
            sim.integrate(t)
        after = [(p.a, p.e, p.inc) for p in sim.particles[1:]]

        self.assertEqual(sim.t, TMAX, msg="%s did not advance time" % integrator)
        for i, (b, a) in enumerate(zip(before, after), start=1):
            for element, x0, x1 in zip(("a", "e", "inc"), b, a):
                with self.subTest(integrator=integrator, particle=i, element=element):
                    self.assertEqual(
                        x0, x1,
                        msg="%s changed particle %d %s from %.17e to %.17e"
                            % (integrator, i, element, x0, x1))

    # -- integrators that damp correctly -----------------------------------
    def test_ias15(self):
        self.check_damped("ias15")

    def test_bs(self):
        self.check_damped("bs")

    def test_whfast(self):
        self.check_damped("whfast")

    def test_saba(self):
        self.check_damped("saba")

    def test_leapfrog(self):
        self.check_damped("leapfrog")

    def test_eos(self):
        self.check_damped("eos")

    # -- integrators that must refuse the operator -------------------------
    def test_mercurius_raises(self):
        self.check_unsupported("mercurius")

    def test_trace_raises(self):
        self.check_unsupported("trace")

    def test_janus_raises(self):
        self.check_unsupported("janus")

    def test_sei_raises(self):
        self.check_unsupported("sei")

    def test_whfast512_raises(self):
        self.check_unsupported("whfast512")

    # -- integrators that do nothing ---------------------------------------
    def test_none_does_nothing(self):
        self.check_noop("none")

if __name__ == '__main__':
    unittest.main()
