#!/bin/bash

#$ -N dihed180
#$ -q tw,pub64,free64
# $ -q free64
#$ -pe openmp 64
#$ -ckpt blcr


# Load the Multicore (OpenMP) NAMD version 2.11:
module load  NAMD-multicore/2.11

# We define $CORES in Grid Engine to make life easier.  The $CORES
# environment variable contains the number of cores allocated
# to this job.

# Replace the input / output files:
charmrun +p $CORES  `which namd2` equil1.inp > out1.log
