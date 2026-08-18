#!/usr/bin/env bash

#SBATCH --time=330:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=embryo_ias15_10-15
#SBATCH --output=embryo_ias15_8.out
#SBATCH --partition=ses



# COMMANDS TO RUN

cd /home/tajer.1/reboundx/examples/embryo_ias15
make
wait
./rebound $SLURM_ARRAY_TASK_ID
wait

