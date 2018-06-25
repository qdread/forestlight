#!/bin/bash --login
#PBS -l walltime=48:00:00,nodes=1:ppn=1,mem=4gb
#PBS -N fitlight
#PBS -j oe
#PBS -m n
#PBS -t 1-3

# Include guild and year as variables with -v flag
# chain is the PBS_ARRAYID

NS=1000
NW=3000

module load GNU/6.2

~/forestlight/stancode/vonb sample algorithm=hmc engine=nuts max_depth=20 num_samples=${NS} num_warmup=${NW} thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_light_${guild}_${year}.r output file=~/forestlight/stanoutput/fit_light_${guild}_${year}_${PBS_ARRAYID}.csv
