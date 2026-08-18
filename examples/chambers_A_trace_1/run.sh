#!/usr/bin/env bash

#SBATCH --time=330:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=cham_org_1
#SBATCH --output=slurm.out
#SBATCH --partition=ses



# COMMANDS TO RUN

cd /home/tajer.1/rebound/examples/cham_org_5
make
wait
./rebound 5
wait
