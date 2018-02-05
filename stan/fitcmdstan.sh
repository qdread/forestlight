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
	~/forestlight/stancode/model_pareto_x_power sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r init=0.1 output file=~/forestlight/stanoutput/fit_paretoxpower_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "weibull" ]; then
	~/forestlight/stancode/model_weibull_x_powerexp sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r init=0.1 output file=~/forestlight/stanoutput/fit_weibullxexp_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "weibullsub" ]; then
	~/forestlight/stancode/model_weibull_x_powerexp sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/ssdump_${guild}_${year}.r init=0.1 output file=~/forestlight/stanoutput/ssfit_weibullxexp_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "paretoexp" ]; then
	~/forestlight/stancode/model_pareto_x_powerexp sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r init=0.1 output file=~/forestlight/stanoutput/fit_paretoxexp_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "weibullpower" ]; then
	~/forestlight/stancode/model_weibull_x_power sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r init=0.1 output file=~/forestlight/stanoutput/fit_weibullxpower_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "paretoexpsub" ]; then
	~/forestlight/stancode/model_pareto_x_powerexp sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/ssdump_${guild}_${year}.r init=0.1 output file=~/forestlight/stanoutput/ssfit_paretoxexp_${guild}_${year}_${PBS_ARRAYID}.csv
fi

if [ "$model" == "weibullpowersub" ]; then
	~/forestlight/stancode/model_weibull_x_power sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/ssdump_${guild}_${year}.r init=0.1 output file=~/forestlight/stanoutput/ssfit_weibullxpower_${guild}_${year}_${PBS_ARRAYID}.csv
fi
