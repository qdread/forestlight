#!/bin/bash --login
#PBS -l walltime=48:00:00,nodes=1:ppn=1,mem=4gb
#PBS -N stanfit
#PBS -j oe
#PBS -m n
#PBS -t 1-3

# Include model as variable with -v flag
# Include guild and year as variables with -v flag
# chain is the PBS_ARRAYID

module load GNU/6.2

if [ "$model" == "pareto" ]; then
	~/forestlight/stancode/model_ppow sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r init=~/forestlight/stanrdump/init_ppow.R output file=~/forestlight/stanoutput/fit_paretoxpower_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "weibull" ]; then
	~/forestlight/stancode/model_wexp sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r init=~/forestlight/stanrdump/init_wpow.R output file=~/forestlight/stanoutput/fit_weibullxexp_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "weibullsub" ]; then
	~/forestlight/stancode/model_wexp sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/ssdump_${guild}_${year}.r init=~/forestlight/stanrdump/init_wexp.R output file=~/forestlight/stanoutput/ssfit_weibullxexp_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "paretoexp" ]; then
	~/forestlight/stancode/model_pexp sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r init=~/forestlight/stanrdump/init_pexp.R output file=~/forestlight/stanoutput/fit_paretoxexp_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "weibullpower" ]; then
	~/forestlight/stancode/model_wpow sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r init=~/forestlight/stanrdump/init_wpow.R output file=~/forestlight/stanoutput/fit_weibullxpower_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "paretoexpsub" ]; then
	~/forestlight/stancode/model_pexp sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/ssdump_${guild}_${year}.r init=~/forestlight/stanrdump/init_pexp.R output file=~/forestlight/stanoutput/ssfit_paretoxexp_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "weibullpowersub" ]; then
	~/forestlight/stancode/model_wpow sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/ssdump_${guild}_${year}.r init=~/forestlight/stanrdump/init_wpow.R output file=~/forestlight/stanoutput/ssfit_weibullxpower_${guild}_${year}_${PBS_ARRAYID}.csv
fi
