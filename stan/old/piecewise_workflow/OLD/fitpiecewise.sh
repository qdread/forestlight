#!/bin/bash --login
#PBS -l walltime=48:00:00,nodes=1:ppn=1,mem=4gb
#PBS -N stanfit
#PBS -j oe
#PBS -m n
#PBS -t 1-3

# Include densitymodel and productionmodel as variables with -v flag
# Include guild and year as variables with -v flag
# Include sampling (NS) and warmup (NW) numbers with -v flag
# Include dumptype (dump for the full data dump, ssdump for the subset dump) with -v flag
# chain is the PBS_ARRAYID

module load GNU/6.2

~/forestlight/stancode/model_d${densitymodel}p${productionmodel} sample algorithm=hmc engine=nuts max_depth=20 num_samples=${NS} num_warmup=${NW} thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/${dumptype}_${guild}_${year}.r output file=~/forestlight/stanoutput/fit_d${densitymodel}p${productionmodel}_${guild}_${year}_${PBS_ARRAYID}.csv

