# Piecewise fits : scaling with light per area as the x variable.
# Converted from qsub to sbatch
# All done with random seeds and proper number of iterations

# D1P1

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=1,productionmodel=1,NS=1000,NW=5000,seed=11 --time=4:00:00 --job-name=d1p1_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=1,productionmodel=1,NS=1000,NW=5000,seed=12 --time=4:00:00 --job-name=d1p1_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=1,productionmodel=1,NS=1000,NW=5000,seed=13 --time=4:00:00 --job-name=d1p1_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=1,productionmodel=1,NS=1000,NW=5000,seed=14 --time=4:00:00 --job-name=d1p1_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=1,productionmodel=1,NS=1000,NW=5000,seed=15 --time=4:00:00 --job-name=d1p1_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=1,productionmodel=1,NS=1000,NW=5000,seed=16 --time=4:00:00 --job-name=d1p1_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=1,productionmodel=1,NS=1000,NW=5000,seed=17 --time=4:00:00 --job-name=d1p1_95_a fitpiecewiselight.sb

# D1P2

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=1,productionmodel=2,NS=1000,NW=5000,seed=21 --time=4:00:00 --job-name=d1p2_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=1,productionmodel=2,NS=1000,NW=5000,seed=22 --time=4:00:00 --job-name=d1p2_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=1,productionmodel=2,NS=1000,NW=5000,seed=23 --time=7:00:00:00 --job-name=d1p2_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=1,productionmodel=2,NS=1000,NW=5000,seed=24 --time=7:00:00:00 --job-name=d1p2_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=1,productionmodel=2,NS=1000,NW=5000,seed=25 --time=4:00:00 --job-name=d1p2_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=1,productionmodel=2,NS=1000,NW=5000,seed=26 --time=4:00:00 --job-name=d1p2_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=1,productionmodel=2,NS=1000,NW=5000,seed=27 --time=7:00:00:00 --job-name=d1p2_95_a fitpiecewiselight.sb

# D2P1

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=31 --time=4:00:00 --job-name=d2p1_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=32 --time=4:00:00 --job-name=d2p1_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=33 --time=2:00:00:00 --job-name=d2p1_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=34 --time=2:00:00:00 --job-name=d2p1_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=35 --time=4:00:00 --job-name=d2p1_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=36 --time=4:00:00 --job-name=d2p1_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=37 --time=2:00:00:00 --job-name=d2p1_95_a fitpiecewiselight.sb

# D2P2

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=41 --time=4:00:00 --job-name=d2p2_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=42 --time=7:00:00:00 --job-name=d2p2_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=43 --time=7:00:00:00 --job-name=d2p2_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=44 --time=3:00:00:00 --job-name=d2p2_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=45 --time=2:00:00:00 --job-name=d2p2_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=46 --time=4:00:00 --job-name=d2p2_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=47 --time=7:00:00:00 --job-name=d2p2_95_a fitpiecewiselight.sb

# D3P1

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=51 --time=2:00:00:00 --job-name=d3p1_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=52 --time=4:00:00 --job-name=d3p1_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=53 --time=7:00:00:00 --job-name=d3p1_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=54 --time=2:00:00:00 --job-name=d3p1_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=55 --time=4:00:00 --job-name=d3p1_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=56 --time=4:00:00 --job-name=d3p1_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=57 --time=7:00:00:00 --job-name=d3p1_95_a fitpiecewiselight.sb

# D3P2

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=61 --time=4:00:00 --job-name=d3p2_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=62 --time=2:00:00:00 --job-name=d3p2_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=63 --time=7:00:00:00 --job-name=d3p2_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=64 --time=7:00:00:00 --job-name=d3p2_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=65 --time=4:00:00:00 --job-name=d3p2_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=66 --time=4:00:00 --job-name=d3p2_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=67 --time=7:00:00:00 --job-name=d3p2_95_a fitpiecewiselight.sb

