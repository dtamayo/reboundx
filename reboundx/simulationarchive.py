import rebound
import reboundx

class Simulationarchive(rebound.Simulationarchive):
    """
    Simulationarchive Class.
    """
    def __new__(cls, filename, rebxfilename, *args, **kwargs):
        return super(Simulationarchive, cls).__new__(cls, filename, *args, **kwargs)

    def __init__(self, filename, rebxfilename, *args, **kwargs):
        """
        Arguments
        ---------
        filename : str
            Filename of the Simulationarchive file to be opened.
        rebxfilename : str
            Filename of the REBOUNDx binary file.
        """
        super(Simulationarchive, self).__init__(filename, *args, **kwargs)
        self.rebxfilename = rebxfilename
        sim, rebx = self[0] # test you can open rebxfilename to warn user if not

    def __getitem__(self, key):
        sim = super(Simulationarchive, self).__getitem__(key)
        rebx = reboundx.Extras(sim, self.rebxfilename)
        return sim, rebx
    
    def getSimulation(self, t, mode='snapshot', keep_unsynchronized=1):
        # same as REBOUND's getSimulation, but adds rebx before integrating to nearby time
        if mode not in ['snapshot', 'close', 'exact']:
            raise AttributeError("Unknown mode.")

        bi, bt = self._getSnapshotIndex(t)
        sim = rebound.Simulation()
        w = c_int(0)
        clibrebound.reb_simulation_init_from_simulationarchive_with_messages(
                byref(sim), byref(self), c_int64(bi), byref(w))

        rebx = reboundx.Extras(sim, self.rebxfilename)

        try:
            safe_mode = sim.integrator.safe_mode
            if safe_mode == 1 or mode == 'exact':
                keep_unsynchronized = 0
            sim.integrator.keep_unsynchronized = keep_unsynchronized
        except AttributeError:
            pass  # not all integrators support keep_unsynchronized

        if mode == 'snapshot':
            sim.synchronize()
        if mode == 'close':
            sim.integrate(t, exact_finish_time=0)
        if mode == 'exact':
            sim.integrate(t, exact_finish_time=1)

        return sim, rebx

    def getBezierPaths(self, origin=None):
        raise AttributeError("getBezierPaths not implemented with REBOUNDx Simulationarchives")
