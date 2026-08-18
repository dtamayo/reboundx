#!/usr/bin/env bash

#SBATCH --time=330:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=embryo_trace_Jul12
#SBATCH --output=embryo_trace_Jul12_%A_%a.out
#SBATCH --partition=ast



# COMMANDS TO RUN

cd /home/tajer.1/reboundx/examples/embryo_TRACE_Jul22
wait
./rebound $SLURM_ARRAY_TASK_ID
wait

