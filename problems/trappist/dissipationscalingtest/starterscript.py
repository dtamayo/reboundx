import numpy as np
from subprocess import call

mags = [0,1.e-3,1.e-2,1.e-1,1]

for mag in mags:
    call("python makestarter.py {0} &".format(mag), shell=True)
