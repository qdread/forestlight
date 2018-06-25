#!/bin/bash --login
#PBS -l walltime=48:00:00,nodes=1:ppn=1,mem=4gb
#PBS -N stanfit
#PBS -j oe
#PBS -m n
#PBS -t 1-3

# Include model as variable with -v flag
# Include guild and year as variables with -v flag
# Include sampling and warmup numbers with -v flag as well.
# chain is the PBS_ARRAYID

module load GNU/6.2

if [ "$model" == "paretopow" ]; then
	~/forestlight/stancode/model_ppow_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=${NS} num_warmup=${NW} thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r output file=~/forestlight/stanoutput/fit_paretoxpower_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "weibullexp" ]; then
	~/forestlight/stancode/model_wexp_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=${NS} num_warmup=${NW} thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r output file=~/forestlight/stanoutput/fit_weibullxexp_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "weibullexpsub" ]; then
	~/forestlight/stancode/model_wexp_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=${NS} num_warmup=${NW} thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/ssdump_${guild}_${year}.r output file=~/forestlight/stanoutput/ssfit_weibullxexp_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "paretoexp" ]; then
	~/forestlight/stancode/model_pexp_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=${NS} num_warmup=${NW} thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r output file=~/forestlight/stanoutput/fit_paretoxexp_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "weibullpow" ]; then
	~/forestlight/stancode/model_wpow_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=${NS} num_warmup=${NW} thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r output file=~/forestlight/stanoutput/fit_weibullxpower_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "paretoexpsub" ]; then
	~/forestlight/stancode/model_pexp_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=${NS} num_warmup=${NW} thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/ssdump_${guild}_${year}.r output file=~/forestlight/stanoutput/ssfit_paretoxexp_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "weibullpowsub" ]; then
	~/forestlight/stancode/model_wpow_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=${NS} num_warmup=${NW} thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/ssdump_${guild}_${year}.r output file=~/forestlight/stanoutput/ssfit_weibullxpower_${guild}_${year}_${PBS_ARRAYID}.csv
fi
