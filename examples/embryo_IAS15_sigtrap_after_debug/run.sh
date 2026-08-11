#!/usr/bin/env bash

#SBATCH --time=330:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=embryo_ias15_sigtrap
#SBATCH --output=embryo_ias15_sigtrap_%A_%a.out



# COMMANDS TO RUN

cd /home/tajer.1/reboundx/examples/embryo_IAS15_sigtrap
wait
./rebound 1
wait

