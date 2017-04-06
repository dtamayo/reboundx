import numpy as np
import sys
from subprocess import call
import time

tmax = 1.e7
Nruns = int(sys.argv[1])

'''
p1 = subprocess.Popen(shlex.split('qstat -a'), stdout=subprocess.PIPE)
p2 = subprocess.Popen(shlex.split('grep dtamayo'), stdin=p1.stdout, stdout=subprocess.PIPE)
p3 = subprocess.Popen(shlex.split('wc -l'), stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
so, se = p3.communicate()

Njobs = int(so.decode())
Nruns = 64-Njobs
'''

with open("runlog.txt", "a+") as f:
    f.seek(0)
    lines = f.readlines()
    if lines:
        lastseed = int(lines[-1])
    else:
        lastseed = -1 

    for i in range(lastseed+1, lastseed+Nruns+1):
        f.write("{0}\n".format(i))
        with open("sunnyvale.sh".format(i), "w") as of:
            of.write("#!/bin/bash -l\n")
            of.write("#PBS -l nodes=1:ppn=1\n")
            of.write("#PBS -q greenq\n")
            of.write("#PBS -r n\n")
            of.write("#PBS -l walltime=168:00:00\n")
            of.write("#PBS -N trappist\n")
            of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
            of.write("cd $PBS_O_WORKDIR\n")
            of.write("module load gcc/4.9.2\n")
            of.write("source /mnt/raid-cita/dtamayo/p3/bin/activate\n")
            of.write("python makeIC.py {0}\n".format(i, tmax))
            of.write("python run.py {0} {1}".format(i, tmax))

        call("chmod u=rwx sunnyvale.sh".format(i), shell=True)
        call("qsub -W x=FLAGS:ADVRES:dtamayo.0 sunnyvale.sh".format(i), shell=True)
