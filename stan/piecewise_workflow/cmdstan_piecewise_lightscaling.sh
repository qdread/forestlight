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
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=1,productionmodel=2,NS=1000,NW=5000,seed=23 --time=7-00:00:00 --job-name=d1p2_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=1,productionmodel=2,NS=1000,NW=5000,seed=24 --time=7-00:00:00 --job-name=d1p2_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=1,productionmodel=2,NS=1000,NW=5000,seed=25 --time=4:00:00 --job-name=d1p2_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=1,productionmodel=2,NS=1000,NW=5000,seed=26 --time=4:00:00 --job-name=d1p2_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=1,productionmodel=2,NS=1000,NW=5000,seed=27 --time=7-00:00:00 --job-name=d1p2_95_a fitpiecewiselight.sb

# D2P1

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=31 --time=4:00:00 --job-name=d2p1_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=32 --time=4:00:00 --job-name=d2p1_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=33 --time=2-00:00:00 --job-name=d2p1_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=34 --time=2-00:00:00 --job-name=d2p1_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=35 --time=4:00:00 --job-name=d2p1_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=36 --time=4:00:00 --job-name=d2p1_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=2,productionmodel=1,NS=1000,NW=5000,seed=37 --time=2-00:00:00 --job-name=d2p1_95_a fitpiecewiselight.sb

# D2P2

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=41 --time=4:00:00 --job-name=d2p2_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=10000,seed=42424 --time=7-00:00:00 --job-name=d2p2_95_2 fitpiecewiselight.sb # rerun with different initial cond. because it did not work.
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=43 --time=7-00:00:00 --job-name=d2p2_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=44 --time=3-00:00:00 --job-name=d2p2_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=45 --time=2-00:00:00 --job-name=d2p2_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=46 --time=4:00:00 --job-name=d2p2_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=5000,seed=47 --time=7-00:00:00 --job-name=d2p2_95_a fitpiecewiselight.sb

# DWP1 (Weibull)

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=w,productionmodel=1,NS=1000,NW=5000,seed=51 --time=2-00:00:00 --job-name=dwp1_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=w,productionmodel=1,NS=1000,NW=5000,seed=52 --time=4:00:00 --job-name=dwp1_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=w,productionmodel=1,NS=1000,NW=5000,seed=53 --time=7-00:00:00 --job-name=dwp1_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=w,productionmodel=1,NS=1000,NW=5000,seed=54 --time=2-00:00:00 --job-name=dwp1_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=w,productionmodel=1,NS=1000,NW=5000,seed=55 --time=4:00:00 --job-name=dwp1_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=w,productionmodel=1,NS=1000,NW=5000,seed=56 --time=4:00:00 --job-name=dwp1_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=w,productionmodel=1,NS=1000,NW=5000,seed=57 --time=7-00:00:00 --job-name=dwp1_95_a fitpiecewiselight.sb

# DWP2 (Weibull)

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=w,productionmodel=2,NS=1000,NW=5000,seed=61 --time=4:00:00 --job-name=dwp2_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=w,productionmodel=2,NS=1000,NW=5000,seed=62 --time=2-00:00:00 --job-name=dwp2_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=w,productionmodel=2,NS=1000,NW=5000,seed=63 --time=7-00:00:00 --job-name=dwp2_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=w,productionmodel=2,NS=1000,NW=5000,seed=64 --time=7-00:00:00 --job-name=dwp2_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=w,productionmodel=2,NS=1000,NW=5000,seed=65 --time=4-00:00:00 --job-name=dwp2_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=w,productionmodel=2,NS=1000,NW=5000,seed=66 --time=4:00:00 --job-name=dwp2_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=w,productionmodel=2,NS=1000,NW=5000,seed=67 --time=7-00:00:00 --job-name=dwp2_95_a fitpiecewiselight.sb

# Final redo of #11 tasks 2 and 3
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=2,productionmodel=2,NS=1000,NW=10000,seed=4224 --time=1-00:00:00 --job-name=d2p2_95_2 --array=2-3 fitpiecewiselight.sb # rerun with different initial cond. because it did not work.

# Extraction job array
sbatch --export=script=extract_ci_piecewiselightfits --mem=16gb --array=1-42 --job-name=pwlightextract jobarray.sb

################ Additional density distributions: 3 piece and lognormal

# D3P1

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=71 --time=2-00:00:00 --job-name=d3p1_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=72 --time=4:00:00 --job-name=d3p1_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=73 --time=7-00:00:00 --job-name=d3p1_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=74 --time=2-00:00:00 --job-name=d3p1_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=75 --time=4:00:00 --job-name=d3p1_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=76 --time=4:00:00 --job-name=d3p1_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=3,productionmodel=1,NS=1000,NW=5000,seed=77 --time=7-00:00:00 --job-name=d3p1_95_a fitpiecewiselight.sb

# D3P2

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=81 --time=4:00:00 --job-name=d3p2_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=82 --time=2-00:00:00 --job-name=d3p2_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=83 --time=7-00:00:00 --job-name=d3p2_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=84 --time=7-00:00:00 --job-name=d3p2_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=85 --time=4-00:00:00 --job-name=d3p2_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=86 --time=4:00:00 --job-name=d3p2_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=3,productionmodel=2,NS=1000,NW=5000,seed=87 --time=7-00:00:00 --job-name=d3p2_95_a fitpiecewiselight.sb

# DlnP1

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=ln,productionmodel=1,NS=1000,NW=5000,seed=91 --time=2-00:00:00 --job-name=dlnp1_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=ln,productionmodel=1,NS=1000,NW=5000,seed=92 --time=4:00:00 --job-name=dlnp1_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=ln,productionmodel=1,NS=1000,NW=5000,seed=93 --time=7-00:00:00 --job-name=dlnp1_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=ln,productionmodel=1,NS=1000,NW=5000,seed=94 --time=2-00:00:00 --job-name=dlnp1_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=ln,productionmodel=1,NS=1000,NW=5000,seed=95 --time=4:00:00 --job-name=dlnp1_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=ln,productionmodel=1,NS=1000,NW=5000,seed=96 --time=4:00:00 --job-name=dlnp1_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=ln,productionmodel=1,NS=1000,NW=5000,seed=97 --time=7-00:00:00 --job-name=dlnp1_95_a fitpiecewiselight.sb

# DlnP2

sbatch --export=dumptype=ssdump,guild=fg1,year=1995,densitymodel=ln,productionmodel=2,NS=1000,NW=5000,seed=101 --time=4:00:00 --job-name=dlnp2_95_1 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg2,year=1995,densitymodel=ln,productionmodel=2,NS=1000,NW=5000,seed=102 --time=2-00:00:00 --job-name=dlnp2_95_2 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg3,year=1995,densitymodel=ln,productionmodel=2,NS=1000,NW=5000,seed=103 --time=7-00:00:00 --job-name=dlnp2_95_3 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg4,year=1995,densitymodel=ln,productionmodel=2,NS=1000,NW=5000,seed=104 --time=7-00:00:00 --job-name=dlnp2_95_4 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=fg5,year=1995,densitymodel=ln,productionmodel=2,NS=1000,NW=5000,seed=105 --time=4-00:00:00 --job-name=dlnp2_95_5 fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=unclassified,year=1995,densitymodel=ln,productionmodel=2,NS=1000,NW=5000,seed=106 --time=4:00:00 --job-name=dlnp2_95_u fitpiecewiselight.sb
sbatch --export=dumptype=ssdump,guild=alltree,year=1995,densitymodel=ln,productionmodel=2,NS=1000,NW=5000,seed=107 --time=7-00:00:00 --job-name=dlnp2_95_a fitpiecewiselight.sb
