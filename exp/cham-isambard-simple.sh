#!/bin/bash
#PBS -q arm
#PBS -l select=1
#PBS -l walltime=00:05:00

RUNCMD="aprun -n 1 -d 64 -j 1  "
#RUNCMD="aprun -n 1 -d 64 -j 1 perf stat "

dim=20000; 
KMP_AFFINITY=disabled  \
    $RUNCMD  \
    $HOME/hicma-dev/chameleon/build/timing/time_dpotrf_tile \
    --threads=64 --M=$dim --N=$dim --K=$dim \
    --P=1 -b 300

    
#chameleon dgemm achieves ~861gflop/s

#$HOME/hicma-dev/chameleon/build/timing/time_dgemm_tile \
