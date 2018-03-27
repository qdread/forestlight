module load GNU/6.2

cd ~/cmdstan-2.15.0
make build
make O=3 ~/forestlight/stancode/model_ppow_withlik
make O=3 ~/forestlight/stancode/model_wpow_withlik
make O=3 ~/forestlight/stancode/model_pexp_withlik
make O=3 ~/forestlight/stancode/model_wexp_withlik

guild="fg1"
year="1995"
PBS_ARRAYID=1

~/forestlight/stancode/model_ppow_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r init=~/forestlight/stanrdump/init_ppow.R output file=~/forestlight/stanoutput/fit_paretoxpower_${guild}_${year}_${PBS_ARRAYID}.csv
~/forestlight/stancode/model_wpow_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r init=~/forestlight/stanrdump/init_wpow.R output file=~/forestlight/stanoutput/fit_weibullxpower_${guild}_${year}_${PBS_ARRAYID}.csv
~/forestlight/stancode/model_pexp_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r init=~/forestlight/stanrdump/init_pexp.R output file=~/forestlight/stanoutput/fit_paretoxexp_${guild}_${year}_${PBS_ARRAYID}.csv
~/forestlight/stancode/model_wexp_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r init=~/forestlight/stanrdump/init_wexp.R output file=~/forestlight/stanoutput/fit_weibullxexp_${guild}_${year}_${PBS_ARRAYID}.csv

~/forestlight/stancode/model_ppow_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r  output file=~/forestlight/stanoutput/fit_paretoxpower_${guild}_${year}_${PBS_ARRAYID}.csv
~/forestlight/stancode/model_wpow_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r  output file=~/forestlight/stanoutput/fit_weibullxpower_${guild}_${year}_${PBS_ARRAYID}.csv
~/forestlight/stancode/model_pexp_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r  output file=~/forestlight/stanoutput/fit_paretoxexp_${guild}_${year}_${PBS_ARRAYID}.csv
~/forestlight/stancode/model_wexp_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r  output file=~/forestlight/stanoutput/fit_weibullxexp_${guild}_${year}_${PBS_ARRAYID}.csv

~/forestlight/stancode/model_ppow_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r  init=0.1 output file=~/forestlight/stanoutput/fit_paretoxpower_${guild}_${year}_${PBS_ARRAYID}.csv
~/forestlight/stancode/model_wpow_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r  init=0.1 output file=~/forestlight/stanoutput/fit_weibullxpower_${guild}_${year}_${PBS_ARRAYID}.csv
~/forestlight/stancode/model_pexp_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r  init=0.1 output file=~/forestlight/stanoutput/fit_paretoxexp_${guild}_${year}_${PBS_ARRAYID}.csv
~/forestlight/stancode/model_wexp_withlik sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dump_${guild}_${year}.r  init=5 output file=~/forestlight/stanoutput/fit_weibullxexp_${guild}_${year}_${PBS_ARRAYID}.csv