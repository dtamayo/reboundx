import numpy as np
import sys
from subprocess import call
import time
import pandas as pd

tmax = 5.e6
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
        with open("sunnyvale.sh".format(i), "w") as of:
            of.write("#!/bin/bash -l\n")
            of.write("#PBS -l nodes=1:ppn=1\n")
            of.write("#PBS -q greenq\n")
            of.write("#PBS -r n\n")
            of.write("#PBS -l walltime=168:00:00\n")
            of.write("#PBS -N trapnodis{0}\n".format(run.name))
            of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
            of.write("cd $PBS_O_WORKDIR\n")
            of.write("module load gcc/4.9.2\n")
            of.write("source /mnt/raid-cita/dtamayo/p3/bin/activate\n")
            of.write("python run.py {0} {1}".format(run['filename'], tmax))

        call("chmod u=rwx sunnyvale.sh", shell=True)
        call("qsub -W x=FLAGS:ADVRES:dtamayo.0 sunnyvale.sh", shell=True)
