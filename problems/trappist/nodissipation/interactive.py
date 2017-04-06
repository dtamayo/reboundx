import numpy as np
import sys
from subprocess import call
import time
import pandas as pd

tmax = 5.e7
Nruns = int(sys.argv[1])

df = pd.read_csv('../data.csv', index_col=0)

with open("data/runlog.txt", "a+") as f:
    f.seek(0)
    lines = f.readlines()
    if lines:
        lastseed = int(lines[-1])
    else:
        lastseed = -1 
    
    for i in range(lastseed+1, lastseed+Nruns+1):
        run = df.iloc[i]
        f.write("{0}\n".format(i))
        print(i)
        call("python run.py {0} {1} &".format(run['filename'], tmax), shell=True)

