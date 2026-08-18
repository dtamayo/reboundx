#!/usr/bin/env bash

#SBATCH --time=330:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=ias15_2_2
#SBATCH --output=ias15_2_2.out
#SBATCH --partition=ast


# COMMANDS TO RUN

cd /home/tajer.1/reboundx/examples/embryo_ias15_2_2
make
wait
./rebound
wait

