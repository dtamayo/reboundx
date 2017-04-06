from subprocess import call
import numpy as np

Nruns=64

for i in range(Nruns):
    with open("sunnyvale.sh", "w") as of:
        of.write("#!/bin/bash -l\n")
        of.write("#PBS -l nodes=1:ppn=1\n")
        of.write("#PBS -q greenq\n")
        of.write("#PBS -r n\n")
        of.write("#PBS -l walltime=48:00:00\n")
        of.write("#PBS -N tides\n")
        of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
        of.write("cd $PBS_O_WORKDIR\n")
        of.write("module load gcc/4.9.2\n")
        of.write("source /mnt/raid-cita/dtamayo/p3/bin/activate\n")
        of.write("python makeIC.py {0}".format(i))

    call("chmod u=rwx sunnyvale.sh", shell=True)
    call("qsub -W x=FLAGS:ADVRES:dtamayo.0 sunnyvale.sh", shell=True)
