# Piecewise fits 
# Converted from qsub to sbatch
# All done with random seeds and proper number of iterations

# D1

sbatch --export=dumptype=dump,guild=fg1,year=1995,model=density1,NS=1000,NW=5000,seed=11 --time=4:00:00 --job-name=d1_1 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg2,year=1995,model=density1,NS=1000,NW=5000,seed=12 --time=4:00:00 --job-name=d1_2 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg3,year=1995,model=density1,NS=1000,NW=5000,seed=13 --time=4:00:00 --job-name=d1_3 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg4,year=1995,model=density1,NS=1000,NW=5000,seed=14 --time=4:00:00 --job-name=d1_4 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg5,year=1995,model=density1,NS=1000,NW=5000,seed=15 --time=4:00:00 --job-name=d1_5 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=unclassified,year=1995,model=density1,NS=1000,NW=5000,seed=16 --time=4:00:00 --job-name=d1_u fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density1,NS=1000,NW=5000,seed=17 --time=4:00:00 --job-name=d1_a fitpiecewiseseparate.sb

# D2

sbatch --export=dumptype=dump,guild=fg1,year=1995,model=density2,NS=1000,NW=5000,seed=211 --time=4:00:00 --job-name=d2_1 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg2,year=1995,model=density2,NS=1000,NW=5000,seed=212 --time=4:00:00 --job-name=d2_2 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg3,year=1995,model=density2,NS=1000,NW=5000,seed=213 --time=7-00:00:00 --job-name=d2_3 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg4,year=1995,model=density2,NS=1000,NW=5000,seed=214 --time=4:00:00 --job-name=d2_4 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg5,year=1995,model=density2,NS=1000,NW=5000,seed=215 --time=4:00:00 --job-name=d2_5 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=unclassified,year=1995,model=density2,NS=1000,NW=5000,seed=216 --time=4:00:00 --job-name=d2_u fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density2,NS=1000,NW=5000,seed=217 --time=7-00:00:00 --job-name=d2_a fitpiecewiseseparate.sb

# D3

sbatch --export=dumptype=dump,guild=fg1,year=1995,model=density3,NS=1000,NW=5000,seed=311 --time=4:00:00 --job-name=d3_1 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg2,year=1995,model=density3,NS=1000,NW=5000,seed=312 --time=4:00:00 --job-name=d3_2 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg3,year=1995,model=density3,NS=1000,NW=5000,seed=313 --time=7-00:00:00 --job-name=d3_3 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg4,year=1995,model=density3,NS=1000,NW=5000,seed=314 --time=4:00:00 --job-name=d3_4 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg5,year=1995,model=density3,NS=1000,NW=5000,seed=315 --time=7-00:00:00 --job-name=d3_5 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=unclassified,year=1995,model=density3,NS=1000,NW=5000,seed=316 --time=4:00:00 --job-name=d3_u fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density3,NS=1000,NW=5000,seed=317 --time=7-00:00:00 --job-name=d3_a fitpiecewiseseparate.sb

# P1

sbatch --export=dumptype=dump,guild=fg1,year=1995,model=production1,NS=1000,NW=5000,seed=411 --time=4:00:00 --job-name=p1_1 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg2,year=1995,model=production1,NS=1000,NW=5000,seed=412 --time=4:00:00 --job-name=p1_2 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg3,year=1995,model=production1,NS=1000,NW=5000,seed=413 --time=4:00:00 --job-name=p1_3 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg4,year=1995,model=production1,NS=1000,NW=5000,seed=414 --time=4:00:00 --job-name=p1_4 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg5,year=1995,model=production1,NS=1000,NW=5000,seed=415 --time=4:00:00 --job-name=p1_5 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=unclassified,year=1995,model=production1,NS=1000,NW=5000,seed=416 --time=4:00:00 --job-name=p1_u fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=production1,NS=1000,NW=5000,seed=417 --time=4:00:00 --job-name=p1_a fitpiecewiseseparate.sb

# P2

sbatch --export=dumptype=dump,guild=fg1,year=1995,model=production2,NS=1000,NW=5000,seed=511 --time=4:00:00 --job-name=p2_1 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg2,year=1995,model=production2,NS=1000,NW=5000,seed=512 --time=4:00:00 --job-name=p2_2 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg3,year=1995,model=production2,NS=1000,NW=5000,seed=513 --time=7-00:00:00 --job-name=p2_3 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg4,year=1995,model=production2,NS=1000,NW=5000,seed=514 --time=7-00:00:00 --job-name=p2_4 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg5,year=1995,model=production2,NS=1000,NW=5000,seed=515 --time=4:00:00 --job-name=p2_5 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=unclassified,year=1995,model=production2,NS=1000,NW=5000,seed=516 --time=4:00:00 --job-name=p2_u fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=production2,NS=1000,NW=5000,seed=517 --time=7-00:00:00 --job-name=p2_a fitpiecewiseseparate.sb

# Additional chains to rerun
# ==========================

# Rerun unfinished jobs (chains that got stuck).
sbatch --export=dumptype=dump,guild=fg5,year=1995,model=density3,NS=1000,NW=5000,seed=3151 --time=7-00:00:00 --job-name=d3_5 --array=1 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density3,NS=1000,NW=5000,seed=666 --time=7-00:00:00 --job-name=d3_a --array=1 fitpiecewiseseparate.sb

# More chains that got stuck (see traceplots)
sbatch --export=dumptype=dump,guild=fg3,year=1995,model=density3,NS=1000,NW=5000,seed=3131 --time=7-00:00:00 --job-name=d3_3 --array=2 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density2,NS=1000,NW=5000,seed=2171 --time=7-00:00:00 --job-name=d2_a fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=fg4,year=1995,model=density3modprior,NS=1000,NW=5000,seed=3141 --time=4:00:00 --job-name=d3_4 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density3,NS=1000,NW=5000,seed=6 --time=7-00:00:00 --job-name=d3_a fitpiecewiseseparate.sb

# This chain just got stuck
sbatch --export=dumptype=dump,guild=fg3,year=1995,model=density3,NS=1000,NW=5000,seed=3132 --time=7-00:00:00 --job-name=d3_3 --array=2 fitpiecewiseseparate.sb

sbatch --export=dumptype=dump,guild=fg4,year=1995,model=density3modprior,NS=1000,NW=10000,seed=6 --time=8:00:00 --job-name=d3_4 fitpiecewiseseparate.sb

sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density2,NS=1000,NW=7500,seed=99 --time=7-00:00:00 --job-name=d2_a fitpiecewiseseparate.sb

sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density3,NS=1000,NW=5000,seed=20 --time=7-00:00:00 --job-name=d3_a --array=3 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density3,NS=1000,NW=5000,seed=61 --time=7-00:00:00 --job-name=d3_a --array=1 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density3,NS=1000,NW=5000,seed=60 --time=7-00:00:00 --job-name=d3_a --array=2 fitpiecewiseseparate.sb

# Runs for June 6th.
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density2,NS=1000,NW=7500,seed=100 --time=7-00:00:00 --job-name=d2_a --array=2 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density3,NS=1000,NW=7500,seed=200 --time=7-00:00:00 --job-name=d3_a --array=1-3 fitpiecewiseseparate.sb

sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density3,NS=1000,NW=7500,seed=208 --time=7-00:00:00 --job-name=d3_a --array=1 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density3,NS=1000,NW=7500,seed=200 --time=7-00:00:00 --job-name=d3_a --array=2 fitpiecewiseseparate.sb
sbatch --export=dumptype=dump,guild=alltree,year=1995,model=density3,NS=1000,NW=7500,seed=204 --time=7-00:00:00 --job-name=d3_a --array=3 fitpiecewiseseparate.sb


# Ensure that the modified prior ones replace the old ones for FG4, D3.
mv fit_density3modprior_fg4_1995_1.csv fit_density3_fg4_1995_1.csv
mv fit_density3modprior_fg4_1995_2.csv fit_density3_fg4_1995_2.csv
mv fit_density3modprior_fg4_1995_3.csv fit_density3_fg4_1995_3.csv

