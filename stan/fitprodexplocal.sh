guilds=("gap" "shade")
years=("1990" "1995" "2000" "2005" "2010" "allyrs")
chains=(1 2 3)

for guild in "${guilds[@]}"; do
	for year in "${years[@]}"; do
		for chain in "${chains[@]}"; do
			~/forestlight/stancode/model_powerlawexp_logtrans sample algorithm=hmc engine=nuts max_depth=20 num_samples=1000 num_warmup=5000 thin=1 adapt delta=0.9 data file=~/forestlight/stanrdump/dat_densprod_${guild}_${year}.r init=0.1 output file=~/forestlight/stanoutput/fit_prod_powerexp_${guild}_${year}_${chain}.csv
		done
	done
done
