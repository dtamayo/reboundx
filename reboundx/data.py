# -*- coding: utf-8 -*-

"""
Initial conditions for standard tests

"""

import math

def add_earths(sim, ei):
    """
    Add a simple test system with 2 ~ Earths, with eccentricities and inclinations (same) specified
    """
    sim.add(m=1.)
    massratio = 3.e-6
    sim.add(m=massratio, a=1., e=ei, inc=ei)
    sim.add(m=massratio, a=3., e=ei, inc=ei)
