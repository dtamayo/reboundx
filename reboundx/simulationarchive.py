import rebound
import reboundx

class SimulationArchive(rebound.SimulationArchive):
    """
    SimulationArchive Class.
    """
    def __init__(self, filename, rebxfilename, *args, **kwargs):
        """
        Arguments
        ---------
        filename : str
            Filename of the SimulationArchive file to be opened.
        rebxfilename : str
            Filename of the REBOUNDx binary file.
        """
        super(SimulationArchive, self).__init__(filename, *args, **kwargs)
        self.rebxfilename = rebxfilename
        sim, rebx = self[0] # test you can open rebxfilename to warn user if not

    def __getitem__(self, key):
        sim = super(SimulationArchive, self).__getitem__(key)
        rebx = reboundx.Extras(sim, self.rebxfilename)
        return sim, rebx

    def getSimulation(self, *args, **kwargs):
        sim = super(SimulationArchive, self).getSimulation(*args, **kwargs)
        rebx = reboundx.Extras(sim, self.rebxfilename)
        return sim, rebx
