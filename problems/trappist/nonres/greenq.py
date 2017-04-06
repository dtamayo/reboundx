import numpy as np
import sys
from subprocess import call
import time
import pandas as pd

Nruns = int(sys.argv[1])
tmax=5.e6
try:
    df = pd.read_csv("../datanonres.csv", index_col=0)
except:
    columns=['filename', 'tinstability']
    df = pd.DataFrame(columns=columns)

lastseed = df.shape[0]

for i in range(lastseed, lastseed+Nruns):
    filename = 'nonresIC{0}.bin'.format(i)
    df.loc[i] = [filename, np.nan]
    with open("sunnyvale.sh".format(i), "w") as of:
        of.write("#!/bin/bash -l\n")
        of.write("#PBS -l nodes=1:ppn=1\n")
        of.write("#PBS -q greenq\n")
        of.write("#PBS -r n\n")
        of.write("#PBS -l walltime=48:00:00\n")
        of.write("#PBS -N trapnonres{0}\n".format(i))
        of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
        of.write("cd $PBS_O_WORKDIR\n")
        of.write("module load gcc/4.9.2\n")
        of.write("source /mnt/raid-cita/dtamayo/p3new/bin/activate\n")
        of.write("python run.py {0} {1} {2}\n".format(i, filename, tmax))

    call("chmod u=rwx sunnyvale.sh", shell=True)
    call("qsub sunnyvale.sh", shell=True)
        
df.to_csv("../datanonres.csv", encoding='ascii')
