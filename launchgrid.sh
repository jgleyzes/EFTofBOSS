#!/bin/bash

#PBS -l nodes=10:ppn=24
#PBS -N COMPUTEGRIDEFT2
#PBS -o /home/pierre/log/
#PBS -e /home/pierre/log/
#PBS -l walltime=1000000000

cd /exports/pierre/EFTofBOSS

mpirun python2 GridChallenge_Mercury.py > /home/pierre/log/grid.log

wait
