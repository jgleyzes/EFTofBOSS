#!/bin/bash
#SBATCH --time=1-16:30:00
#SBATCH --partition=iric,hns,normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5GB

ml load python/2.7.13

python MCMC_newwidow_marg.py $1 $2 $3 $4 $5 $6 $7 $8

