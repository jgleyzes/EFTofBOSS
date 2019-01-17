#!/bin/bash

#PBS -l nodes=10:ppn=24
#PBS -N COMPUTEFSIGMA8
#PBS -o /home/pierre/log/
#PBS -e /home/pierre/log/
#PBS -l walltime=1000000000

cd /exports/pierre/EFTofBOSS

mpirun python2 fsigma8.py > /home/pierre/log/fs8.log

wait
