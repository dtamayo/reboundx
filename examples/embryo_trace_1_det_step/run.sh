#!/usr/bin/env bash

#SBATCH --time=330:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=trace_1
#SBATCH --output=trace_1.out
#SBATCH --partition=ast


# COMMANDS TO RUN

cd /home/tajer.1/reboundx/examples/embryo_trace_1_det_step
wait
make
wait
./rebound 1
wait

