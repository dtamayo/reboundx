#!/usr/bin/env bash

#SBATCH --time=330:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=embryo_trace_0_20
#SBATCH --output=embryo_trace.out



# COMMANDS TO RUN

cd /home/tajer.1/reboundx/examples/embryo_trace
make
wait
./rebound $SLURM_ARRAY_TASK_ID
wait

