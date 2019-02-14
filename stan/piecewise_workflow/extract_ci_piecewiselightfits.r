# 14 Feb 2019: Extract info from piecewise fits scaled by light received.
# For each CMDstan model fit, get the following
# Credible intervals of parameters
# Credible and prediction intervals to make graphs
# Information criteria
# Fitted slopes in log space
# Bayesian R-squared

# New version created 10 Sep: only 1995 for piecewise.
# Edit 10 Aug: specify subset for alltree and for fg 3
# Edit 12 April: Do this in parallel because it is too slow
task <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

source('~/forestlight/stancode/extraction_functions_piecewise.r')

library(purrr)
library(dplyr)
library(rstan, lib.loc = '/mnt/home/qdr/R/x86_64-pc-linux-gnu-library/3.5') # Ensure newest version is loaded.

mod_df <- expand.grid(dens_model = c('1','2','w'),
                      prod_model = c('1','2'),
                      fg = c('fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'alltree', 'unclassified'),
                      year = 1995, 
                      stringsAsFactors = FALSE)

min_n <- read.csv('~/forestlight/stanrdump/min_n_lighttrees.csv', stringsAsFactors = FALSE)

total_prod <- read.csv('~/forestlight/stanrdump/production_total_lighttrees.csv', stringsAsFactors = FALSE)

mod_df <- mod_df %>%
  left_join(total_prod) %>%
  rename(total_production = production) %>%
  left_join(min_n)

# Predicted values.
dbh_pred <- exp(seq(log(1.1), log(412), length.out = 101))

use_subset <- (mod_df$fg[task] %in% c('fg3', 'alltree'))

fit_info <- extract_all_fit(dens_model = mod_df$dens_model[task],
                            prod_model = mod_df$prod_model[task],
                            fg = mod_df$fg[task],
                            year = mod_df$year[task],
                            xmin = mod_df$xmin[task],
                            n = mod_df$n[task],
                            total_production = mod_df$total_production[task],
							use_subset = use_subset,
							fp = '~/forestlight/stanoutput/lightfits',
							fpdump = '~/forestlight/stanrdump',
							fitprefix = 'fit_d',
							dumpprefix = 'ssdump_lightscaling_',
							LL = 1.08,
							UL = 412.2
							)

save(fit_info, file = paste0('~/forestlight/stanoutput/lightfitinfo/pw_info_',task,'.r'))
