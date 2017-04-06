from subprocess import call
import numpy as np

filename = "4res" # filename of start file in data/, without .bin at end
tmax = 5.e6
Nruns = 32 
taue0min = 1.e2
taue0max = 1.e8
taue0s = np.logspace(8, 4, Nruns)

for taue0 in taue0s:
    with open("sunnyvale.sh", "w") as of:
        of.write("#!/bin/bash\n")
        of.write("#PBS -l nodes=1:ppn=1\n")
        of.write("#PBS -q greenq\n")
        of.write("#PBS -r n\n")
        of.write("#PBS -l walltime=24:00:00\n")
        of.write("#PBS -N tides\n")
        of.write("# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n")
        of.write("cd $PBS_O_WORKDIR\n")
        of.write("module load gcc/4.9.2\n")
        of.write("source /mnt/raid-cita/dtamayo/p3/bin/activate\n")
        of.write("python run.py {0} {1} {2}".format(filename, taue0, tmax))

    call("chmod u=rwx sunnyvale.sh", shell=True)
    call("qsub -W x=FLAGS:ADVRES:dtamayo.0 sunnyvale.sh", shell=True)
