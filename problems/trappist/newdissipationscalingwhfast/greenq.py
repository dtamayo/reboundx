import numpy as np
import sys
from subprocess import call
import time

taues = np.logspace(2,3,4)
taues = taues[:-1]
mags = [1.e-3,1.e-1,1]

for mag in mags:
    for i,taue in enumerate(taues):
        with open("sunnyvale.sh".format(i), "w") as of:
            of.write("#!/bin/bash -l\n")
            of.write("#PBS -l nodes=1:ppn=1\n")
            of.write("#PBS -q greenq\n")
            of.write("#PBS -r n\n")
            of.write("#PBS -l walltime=168:00:00\n")
            of.write("#PBS -N trappisttaue{0}\n".format(i))
            of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
            of.write("cd $PBS_O_WORKDIR\n")
            of.write("module load gcc/4.9.2\n")
            of.write("source /mnt/raid-cita/dtamayo/p3/bin/activate\n")
            of.write("python run.py {0} {1}".format(taue, mag))

        call("chmod u=rwx sunnyvale.sh".format(i), shell=True)
        call("qsub -W x=FLAGS:ADVRES:dtamayo.0 sunnyvale.sh".format(i), shell=True)
