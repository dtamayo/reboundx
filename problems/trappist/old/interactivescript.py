from subprocess import call
import numpy as np

filename = "6res" # filename of start file in data/, without .bin at end
tmax = 5.e6
Nruns = 13
taue0s = np.logspace(8, 2, Nruns)

for taue0 in taue0s:
    call("python run.py {0} {1} {2} &".format(filename, taue0, tmax), shell=True)
