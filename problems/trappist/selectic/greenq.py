import numpy as np
import sys
from subprocess import call
import time

Nic = 1000
Ncores = 50 

ids = np.arange(0,Nic,int(Nic/(Ncores-1)))
ids = np.append(ids,Nic)

for i in range(ids.shape[0]-1):
    with open("sunnyvale.sh".format(i), "w") as of:
        of.write("#!/bin/bash -l\n")
        of.write("#PBS -l nodes=1:ppn=1\n")
        of.write("#PBS -q greenq\n")
        of.write("#PBS -r n\n")
        of.write("#PBS -l walltime=168:00:00\n")
        of.write("#PBS -N trap\n".format(i))
        of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
        of.write("cd $PBS_O_WORKDIR\n")
        of.write("module load gcc/4.9.2\n")
        of.write("source /mnt/raid-cita/dtamayo/p3/bin/activate\n")
        of.write("python run.py {0} {1} {2}".format(i,ids[i],ids[i+1]))

    call("chmod u=rwx sunnyvale.sh".format(i), shell=True)
    call("qsub -W x=FLAGS:ADVRES:dtamayo.0 sunnyvale.sh".format(i), shell=True)
