#!/bin/bash --login
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=4gb
#PBS -N prod_powerlaw
#PBS -t 1-36
#PBS -j oe
#PBS -m n

module load GNU/6.2

guildnames=("gap" "shade")
yearnames=("1990" "1995" "2000" "2005" "2010" "allyrs")
chainnames=(1 2 3)

guildnos=(0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1)
yearnos=(0 0 1 1 2 2 3 3 4 4 5 5 0 0 1 1 2 2 3 3 4 4 5 5 0 0 1 1 2 2 3 3 4 4 5 5)
chainnos=(0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2)

n=$((PBS_ARRAYID - 1))

guild=${guildnames[guildnos[n]]}
year=${yearnames[yearnos[n]]}
chain=${chainnames[chainnos[n]]}

~/forestlight/stancode/model_powerlaw_logtrans sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dat_densprod_${guild}_${year}.r init=0.1 output file=~/forestlight/stanoutput/fit_prod_powerlaw_${guild}_${year}_${chain}.csv

